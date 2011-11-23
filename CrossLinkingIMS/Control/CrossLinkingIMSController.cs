﻿using System.Collections.Generic;
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
using ProteinDigestionSimulator;

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
		/// <param name="massToleranceBase">Mass tolerance of instrument, in ppm.</param>
		/// <param name="fastAFile">The FileInfo object for the FASTA file containg all protein sequences you want to search.</param>
		/// <param name="featureFile">The FileInfo object for the LC-IMS-MS features file, created by the LC-IMS-MS Feature Finder. (email Kevin.Crowell@pnnl.gov for more info)</param>
		/// <param name="peaksFile">The FileInfo object for the Isotopic Peaks file, created by DeconTools. (email Gordon.Slysz@pnnl.gov for more info)</param>
		/// <returns>An enumerable of CrossLinkResult objects.</returns>
		public static IEnumerable<CrossLinkResult> Execute(double massToleranceBase, FileInfo fastAFile, FileInfo featureFile, FileInfo peaksFile)
		{
			// Read in FASTA File
			FastAParser fastAParser = new FastAParser(fastAFile.FullName);
			IEnumerable<ISequence> sequenceEnumerable = fastAParser.Parse();

			// Read in LC-IMS-MS Features
			List<LcImsMsFeature> featureList = LcImsMsFeatureReader.ReadFile(featureFile);

			// Read in Isotopic Peaks (not Isotopic Profile)
			List<IsotopicPeak> peakEnumerable = IsotopicPeakReader.ReadFile(peaksFile);

			// Now call the executor that expects the opbjects instead of the file locations
			return Execute(massToleranceBase, sequenceEnumerable, featureList, peakEnumerable);
		}

		/// <summary>
		/// Executes the cross-link search for LC-IMS-TOF data.
		/// </summary>
		/// <param name="massToleranceBase">Mass tolerance of instrument, in ppm.</param>
		/// <param name="proteinSequenceEnumerable">IEnumerable of protein sequences, as a .NET Bio ISequence object.</param>
		/// <param name="featureList">List of LC-IMS-MS Features, as LcImsMsFeature.</param>
		/// <param name="peakList">List of Isotopic Peaks, as IsotopicPeak.</param>
		/// <returns>An enumerable of CrossLinkResult objects.</returns>
		public static IEnumerable<CrossLinkResult> Execute(double massToleranceBase, IEnumerable<ISequence> proteinSequenceEnumerable, List<LcImsMsFeature> featureList, List<IsotopicPeak> peakList)
		{
			// Used for finding Isotopic Profiles in the data
			BasicTFF msFeatureFinder = new BasicTFF();

			List<CrossLink> crossLinkList = new List<CrossLink>();

			// Create CrossLink objects from all of the protein sequences
			foreach (ISequence proteinSequence in proteinSequenceEnumerable)
			{
				string proteinSequenceString = new string(proteinSequence.Select((a => (char)a)).ToArray());
				string proteinId = proteinSequence.ID;

				// Get a List of Peptides from the Protein Sequence
				IEnumerable<clsInSilicoDigest.PeptideInfoClass> peptideList = SequenceUtil.DigestProtein(proteinSequenceString);

				// Find all possible cross links from the peptide list
				IEnumerable<CrossLink> crossLinkEnumerable = CrossLinkUtil.GenerateTheoreticalCrossLinks(peptideList, proteinSequenceString, proteinId);
				crossLinkList.AddRange(crossLinkEnumerable);
			}

			// Sort the CrossLinks by mass so that the results are ordered in a friendly way
			IEnumerable<CrossLink> orderedCrossLinkEnumerable = crossLinkList.OrderBy(x => x.Mass);

			// Sort Feature by mass so we can use binary search
			featureList = featureList.OrderBy(x => x.MassMonoisotopic).ToList();

			// Set up a Feature Comparer to use for binary search later on
			AnonymousComparer<LcImsMsFeature> featureComparer = new AnonymousComparer<LcImsMsFeature>((x, y) => x.MassMonoisotopic.CompareTo(y.MassMonoisotopic));

			// Sort the Isotopic Peaks by LC Scan, IMS Scan, and m/z to set them up for binary search later on
			peakList = peakList.OrderBy(x => x.ScanLc).ThenBy(x => x.ScanIms).ThenBy(x => x.Mz).ToList();

			List<CrossLinkResult> crossLinkResultList = new List<CrossLinkResult>();

			// Search the data for the existence of cross-links
			foreach (CrossLink crossLink in orderedCrossLinkEnumerable)
			{
				// Calculate mass tolerance to use for binary search
				double massTolerance = massToleranceBase * crossLink.Mass / GeneralConstants.PPM_DIVISOR;

				LcImsMsFeature lowFeature = new LcImsMsFeature { MassMonoisotopic = crossLink.Mass - massTolerance };
				LcImsMsFeature highFeature = new LcImsMsFeature { MassMonoisotopic = crossLink.Mass + massTolerance };

				int lowFeaturePosition = featureList.BinarySearch(lowFeature, featureComparer);
				int highFeaturePosition = featureList.BinarySearch(highFeature, featureComparer);

				lowFeaturePosition = lowFeaturePosition < 0 ? ~lowFeaturePosition : lowFeaturePosition;
				highFeaturePosition = highFeaturePosition < 0 ? ~highFeaturePosition : highFeaturePosition;

				// Iterate over all LC-IMS-MS Features that match the Unmodified cross-link mass
				for (int i = lowFeaturePosition; i < highFeaturePosition; i++)
				{
					LcImsMsFeature feature = featureList[i];

					// Search for a mass shift in each of the LC Scans the unmodified cross-link mass was found
					for (int currentScanLc = feature.ScanLcStart; currentScanLc <= feature.ScanLcEnd; currentScanLc++)
					{
						CrossLinkResult crossLinkResult = new CrossLinkResult(crossLink, feature, currentScanLc);

						List<IPeak> candidatePeaks = PeakUtil.FindCandidatePeaks(peakList, feature.MzMonoisotopic, currentScanLc, feature.ScanImsRep);
						List<double> massShiftList = crossLink.MassShiftList;
						List<double> shiftedMassList = new List<double>();

						// Calculate the shifted mass values that we want to search for
						switch (massShiftList.Count)
						{
							case 1:
								{
									double firstNewMass = feature.MassMonoisotopic + massShiftList[0];
									shiftedMassList.Add(firstNewMass);
								}
								break;
							case 2:
								{
									double firstNewMass = feature.MassMonoisotopic + massShiftList[0];
									double secondNewMass = feature.MassMonoisotopic + massShiftList[1];
									double thirdNewMass = feature.MassMonoisotopic + massShiftList[0] + massShiftList[1];

									shiftedMassList.Add(firstNewMass);
									shiftedMassList.Add(secondNewMass);
									shiftedMassList.Add(thirdNewMass);
								}
								break;
						}

						// Search for shifted mass values in Isotopic Peaks
						foreach (double shiftedMass in shiftedMassList)
						{
							double shiftedMz = (shiftedMass / feature.ChargeState) + GeneralConstants.MASS_OF_PROTON;

							// Create theoretical Isotopic Peaks that will later form a theoretical Isotopic Profile
							List<MSPeak> theoreticalPeakList = new List<MSPeak> { new MSPeak { XValue = shiftedMz, Height = 1 } };
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
							IsotopicProfile isotopicProfile = new IsotopicProfile
							{
								MonoIsotopicMass = shiftedMass,
								MonoPeakMZ = shiftedMz,
								ChargeState = feature.ChargeState,
								Peaklist = sortPeaksQuery.ToList()
							};

							// Search for the theoretical Isotopic Profile
							IsotopicProfile foundProfile = msFeatureFinder.FindMSFeature(candidatePeaks, isotopicProfile, massToleranceBase, false);

							/*
							 * It is possible that the set mono pass of the previous theoretical distribution was the right-most peak of the actual distribution
							 * If so, we should be able to shift the theoretical distribution over to the left and find the actual distribution
							 */
							if (foundProfile == null)
							{
								foreach (MSPeak msPeak in sortPeaksQuery)
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
							bool didFindProfile = foundProfile != null;
							crossLinkResult.MassShiftResults.KvpList.Add(new KeyValuePair<double, bool>(shiftedMass, didFindProfile));
						}

						crossLinkResultList.Add(crossLinkResult);
					}
				}
			}

			return crossLinkResultList;
		}
	}
}
