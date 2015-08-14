using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Bio;
using Bio.IO.FastA;
using CrossLinkingIMS.Constants;
using CrossLinkingIMS.Data;
using CrossLinkingIMS.IO;
using CrossLinkingIMS.Util;
using DeconTools.Backend.Core;
using DeconTools.Backend.ProcessingTasks.TargetedFeatureFinders;

namespace CrossLinkingIMS.Control
{
    /// <summary>
    /// A class for various ways of executing the cross-linking IMS code.
    /// </summary>
    public class CrossLinkingImsController
    {
        /// <summary>
        /// Executes the cross-link search for LC-IMS-TOF data.
        /// </summary>
        /// <param name="settings">Settings object to control parameters for cross-linking.</param>
        /// <param name="fastAFile">The FileInfo object for the FASTA file containg all protein sequences you want to search.</param>
        /// <param name="featureFile">The FileInfo object for the LC-IMS-MS features file, created by the LC-IMS-MS Feature Finder. (email Kevin.Crowell@pnnl.gov for more info)</param>
        /// <param name="peaksFile">The FileInfo object for the Isotopic Peaks file, created by DeconTools. (email Gordon.Slysz@pnnl.gov for more info)</param>
        /// <returns>An enumerable of CrossLinkResult objects.</returns>
        public static IList<CrossLinkResult> Execute(CrossLinkSettings settings, FileInfo fastAFile, FileInfo featureFile, FileInfo peaksFile)
        {

            IEnumerable<ISequence> sequenceEnumerable;
            List<LcImsMsFeature> featureList;
            List<IsotopicPeak> peakEnumerable;
 
            try
            {

                // Read in FASTA File
                var fastAParser = new FastAParser(fastAFile.FullName);
                sequenceEnumerable = fastAParser.Parse();
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error reading the FASTA file: " + ex.Message);
                throw;
            }

            try
            {
                // Read in LC-IMS-MS Features
                featureList = LcImsMsFeatureReader.ReadFile(featureFile);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error reading the LCMSFeatures file: " + ex.Message);
                throw;
            }

            try
            {
                // Read in Isotopic Peaks (not Isotopic Profile)
                peakEnumerable = IsotopicPeakReader.ReadFile(peaksFile);
            }
            catch (Exception ex)
            {
                Console.WriteLine("Error reading the Isotopic Peaks file: " + ex.Message);
                throw;
            }
            
            // Now call the executor that expects the opbjects instead of the file locations
            return Execute(settings, sequenceEnumerable, featureList, peakEnumerable);
        }

        /// <summary>
        /// Executes the cross-link search for LC-IMS-TOF data.
        /// </summary>
        /// <param name="settings">Settings object to control parameters for cross-linking.</param>
        /// <param name="proteinSequenceEnumerable">IEnumerable of protein sequences, as a .NET Bio ISequence object.</param>
        /// <param name="featureList">List of LC-IMS-MS Features, as LcImsMsFeature.</param>
        /// <param name="peakList">List of Isotopic Peaks, as IsotopicPeak.</param>
        /// <returns>An enumerable of CrossLinkResult objects.</returns>
        public static IList<CrossLinkResult> Execute(
            CrossLinkSettings settings, 
            IEnumerable<ISequence> proteinSequenceEnumerable, 
            List<LcImsMsFeature> featureList, 
            List<IsotopicPeak> peakList)
        {
            var massToleranceBase = settings.MassTolerance;
            var maxMissedCleavages = settings.MaxMissedCleavages;
            var digestionRule = settings.TrypticType;

            CrossLinkUtil.StaticDeltaMass = settings.StaticDeltaMass;
            CrossLinkUtil.UseC13 = settings.UseC13;
            CrossLinkUtil.UseN15 = settings.UseN15;

            // Used for finding Isotopic Profiles in the data
            var msFeatureFinder = new BasicTFF();

            var crossLinkList = new List<CrossLink>();
            var lastProgress = DateTime.UtcNow;
            var proteinsProcessed = 0;

            // Create CrossLink objects from all of the protein sequences
            foreach (var proteinSequence in proteinSequenceEnumerable)
            {
                var proteinSequenceString = new string(proteinSequence.Select((a => (char)a)).ToArray());
                var proteinId = proteinSequence.ID;

                // Get a List of Peptides from the Protein Sequence
                var peptideList = SequenceUtil.DigestProtein(proteinSequenceString, digestionRule, maxMissedCleavages);

                // Find all possible cross links from the peptide list
                var crossLinkEnumerable = CrossLinkUtil.GenerateTheoreticalCrossLinks(peptideList, proteinSequenceString, proteinId);
                crossLinkList.AddRange(crossLinkEnumerable);

                proteinsProcessed++;
                if (DateTime.UtcNow.Subtract(lastProgress).TotalSeconds >= 15)
                {
                    lastProgress = DateTime.UtcNow;
                    Console.WriteLine("Creating cross linked peptide list; " + proteinsProcessed + " proteins processed");
                }

            }

            Console.WriteLine("Sorting cross-linked peptides");

            // Sort the CrossLinks by mass so that the results are ordered in a friendly way
            IEnumerable<CrossLink> orderedCrossLinkEnumerable = crossLinkList.OrderBy(x => x.Mass);

            // Sort Feature by mass so we can use binary search
            featureList = featureList.OrderBy(x => x.MassMonoisotopic).ToList();

            // Set up a Feature Comparer and Peak Comparer to use for binary search later on
            var featureComparer = new AnonymousComparer<LcImsMsFeature>((x, y) => x.MassMonoisotopic.CompareTo(y.MassMonoisotopic));
            var peakComparer = new AnonymousComparer<IsotopicPeak>((x, y) => x.ScanLc != y.ScanLc ? x.ScanLc.CompareTo(y.ScanLc) : x.ScanIms != y.ScanIms ? x.ScanIms.CompareTo(y.ScanIms) : x.Mz.CompareTo(y.Mz));

            // Sort the Isotopic Peaks by LC Scan, IMS Scan, and m/z to set them up for binary search later on
            peakList.Sort(peakComparer);

            var crossLinkResultList = new List<CrossLinkResult>();
            var totalCandidatePeptides = crossLinkList.Count;

            Console.WriteLine("Searching isotopic data vs. " + totalCandidatePeptides.ToString("#,##0") + " candidate cross-linked peptides");
            lastProgress = DateTime.UtcNow;
            var crosslinkCandidatesProcessed = 0;
 
            // Search the data for the existence of cross-links
            foreach (var crossLink in orderedCrossLinkEnumerable)
            {
                // Calculate mass tolerance to use for binary search
                var massTolerance = massToleranceBase * crossLink.Mass / GeneralConstants.PPM_DIVISOR;

                var lowFeature = new LcImsMsFeature { MassMonoisotopic = crossLink.Mass - massTolerance };
                var highFeature = new LcImsMsFeature { MassMonoisotopic = crossLink.Mass + massTolerance };

                var lowFeaturePosition = featureList.BinarySearch(lowFeature, featureComparer);
                var highFeaturePosition = featureList.BinarySearch(highFeature, featureComparer);

                lowFeaturePosition = lowFeaturePosition < 0 ? ~lowFeaturePosition : lowFeaturePosition;
                highFeaturePosition = highFeaturePosition < 0 ? ~highFeaturePosition : highFeaturePosition;

                // Iterate over all LC-IMS-MS Features that match the Unmodified cross-link mass
                for (var i = lowFeaturePosition; i < highFeaturePosition; i++)
                {
                    var feature = featureList[i];

                    // Search for a mass shift in each of the LC Scans the unmodified cross-link mass was found
                    for (var currentScanLc = feature.ScanLcStart; currentScanLc <= feature.ScanLcEnd; currentScanLc++)
                    {
                        var crossLinkResult = new CrossLinkResult(crossLink, feature, currentScanLc);

                        var candidatePeaks = PeakUtil.FindCandidatePeaks(peakList, feature.MzMonoisotopic, currentScanLc, feature.ScanImsRep);
                        var massShiftList = crossLink.MassShiftList;
                        var shiftedMassList = new List<double>();

                        // Calculate the shifted mass values that we want to search for
                        switch (massShiftList.Count)
                        {
                            case 1:
                                {
                                    var firstNewMass = feature.MassMonoisotopic + massShiftList[0];
                                    shiftedMassList.Add(firstNewMass);
                                }
                                break;
                            case 2:
                                {
                                    var firstNewMass = feature.MassMonoisotopic + massShiftList[0];
                                    var secondNewMass = feature.MassMonoisotopic + massShiftList[1];
                                    var thirdNewMass = feature.MassMonoisotopic + massShiftList[0] + massShiftList[1];

                                    shiftedMassList.Add(firstNewMass);
                                    shiftedMassList.Add(secondNewMass);
                                    shiftedMassList.Add(thirdNewMass);
                                }
                                break;
                        }

                        // Search for shifted mass values in Isotopic Peaks
                        foreach (var shiftedMass in shiftedMassList)
                        {
                            var shiftedMz = (shiftedMass / feature.ChargeState) + GeneralConstants.MASS_OF_PROTON;

                            // Create theoretical Isotopic Peaks that will later form a theoretical Isotopic Profile
                            var theoreticalPeakList = new List<MSPeak> { new MSPeak { XValue = shiftedMz, Height = 1 } };
                            for (double k = 1; k < 4; k++)
                            {
                                theoreticalPeakList.Add(new MSPeak { XValue = shiftedMz + (k * 1.003 / feature.ChargeState), Height = (float)(1.0 - (k / 4)) });
                                theoreticalPeakList.Add(new MSPeak { XValue = shiftedMz - (k * 1.003 / feature.ChargeState), Height = (float)(1.0 - (k / 4)) });
                            }

                            // Sort peaks by m/z
                            var sortPeaksQuery = from peak in theoreticalPeakList
                                                 orderby peak.XValue
                                                 select peak;

                            // Create a theoretical Isotopic Profile for DeconTools to search for
                            var isotopicProfile = new IsotopicProfile
                            {
                                MonoIsotopicMass = shiftedMass,
                                MonoPeakMZ = shiftedMz,
                                ChargeState = feature.ChargeState,
                                Peaklist = sortPeaksQuery.ToList()
                            };

                            // Search for the theoretical Isotopic Profile
                            var foundProfile = msFeatureFinder.FindMSFeature(candidatePeaks, isotopicProfile, massToleranceBase, false);

                            /*
                             * It is possible that the set mono pass of the previous theoretical distribution was the right-most peak of the actual distribution
                             * If so, we should be able to shift the theoretical distribution over to the left and find the actual distribution
                             */
                            if (foundProfile == null)
                            {
                                foreach (var msPeak in sortPeaksQuery)
                                {
                                    msPeak.XValue -= (1.003 / feature.ChargeState);
                                }

                                isotopicProfile = new IsotopicProfile
                                {
                                    MonoIsotopicMass = shiftedMass - 1.003,
                                    MonoPeakMZ = shiftedMz - (1.003 / feature.ChargeState),
                                    ChargeState = feature.ChargeState,
                                    Peaklist = sortPeaksQuery.ToList()
                                };

                                foundProfile = msFeatureFinder.FindMSFeature(candidatePeaks, isotopicProfile, massToleranceBase, false);
                            }

                            // Add to results, even if we did not find it.
                            var didFindProfile = foundProfile != null;
                            crossLinkResult.MassShiftResults.KvpList.Add(new KeyValuePair<double, bool>(shiftedMass, didFindProfile));
                        }

                        crossLinkResultList.Add(crossLinkResult);
                    }
                }

                crosslinkCandidatesProcessed++;
                if (DateTime.UtcNow.Subtract(lastProgress).TotalSeconds >= 10)
                {
                    lastProgress = DateTime.UtcNow;
                    var percentComplete = crosslinkCandidatesProcessed / (double)totalCandidatePeptides * 100;

                    Console.WriteLine("Searching isotopic data; " + percentComplete.ToString("0.0") + "% complete");
                }

            }

            return crossLinkResultList;
        }
    }
}
