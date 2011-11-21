using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
using System.Text;
using CrossLinkingIMS.Data;
using CrossLinkingIMS.IO;
using CrossLinkingIMS.Util;
using DeconTools.Backend.Core;
using DeconTools.Backend.DTO;
using DeconTools.Backend.Data;
using DeconTools.Backend.ProcessingTasks.TargetedFeatureFinders;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS
{
	public class ConsoleApplication
	{
		private static double MASS_TOLERANCE_BASE = 20;
		private const int PPM_DIVISOR = 1000000;
		private const double MASS_OF_PROTON = 1.00727649;

		static void Main(string[] args)
		{
			BasicTFF msFeatureFinder = new BasicTFF();

			// Hard-coded Protein Sequence
			const string proteinSequence = "AEQVSKQEISHFKLVKVGTINVSQSGGQISSPSDLREKLSELADAKGGKYYHIIAAREHGPNFEAVAEVYNDATKLEHHHHHH";

			// Get a List of Peptides from the Protein Sequence
			IEnumerable<clsInSilicoDigest.PeptideInfoClass> peptideList = SequenceUtil.DigestProtein(proteinSequence);

			// Find all possible cross links from the peptide list
			IEnumerable<CrossLink> crossLinkEnumerable = SequenceUtil.GenerateTheoreticalCrossLinks(peptideList, proteinSequence);
			IEnumerable<CrossLink> orderedCrossLinkEnumerable = crossLinkEnumerable.OrderBy(o => o.Mass);

			// Read in LC-IMS-MS Features
			FileInfo featureFile = new FileInfo("SrfN_NC50_IMS_7Sep11_Roc_11-05-30_LCMSFeatures.txt");
			List<LcImsMsFeature> featureList = LcImsMsFeatureReader.ReadFile(featureFile);

			var sortFeatureListQuery = from feature in featureList
			                           orderby feature.MassMonoisotopic
			                           select feature;

			featureList = sortFeatureListQuery.ToList();

			// Set up a Feature Comparer to use for binary search later on
			AnonymousComparer<LcImsMsFeature> featureComparer = new AnonymousComparer<LcImsMsFeature>((x, y) => x.MassMonoisotopic.CompareTo(y.MassMonoisotopic));

			// Read in Isotopic Peaks (not Isotopic Profile)
			FileInfo peaksFile = new FileInfo("SrfN_NC50_IMS_7Sep11_Roc_11-05-30_peaks.txt");
			BackgroundWorker backgroundWorker = new BackgroundWorker {WorkerReportsProgress = true, WorkerSupportsCancellation = true};
			PeakImporterFromText peakImporter = new PeakImporterFromText(peaksFile.FullName, backgroundWorker);
			List<IPeak> iPeakList = new List<IPeak>();
			peakImporter.ImportUIMFPeaks(iPeakList);

			List<MSPeakResult> peakList = iPeakList.Select(i => (MSPeakResult)i).ToList();

			// Sort the Isotopic Peaks by LC Scan, IMS Scan, and m/z to set them up for binary search later on
			var sortPeakListQuery = from peak in peakList
									orderby peak.Frame_num, peak.Scan_num, peak.XValue
									select peak;

			peakList = sortPeakListQuery.ToList();

			List<CrossLinkResult> crossLinkResultList = new List<CrossLinkResult>();

			// Search the data for the existence of cross-links
			foreach (CrossLink crossLink in orderedCrossLinkEnumerable)
			{
				// Calculate mass tolerance to use for binary search
				double massTolerance = MASS_TOLERANCE_BASE * crossLink.Mass / PPM_DIVISOR;

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

						List<IPeak> candidatePeaks = FindCandidatePeaks(peakList, feature.MzMonoisotopic, currentScanLc, feature.ScanImsRep);
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
							double shiftedMz = (shiftedMass/feature.ChargeState) + MASS_OF_PROTON;

							// Create theoretical Isotopic Peaks that will later form a theoretical Isotopic Profile
							List<MSPeak> theoreticalPeakList = new List<MSPeak> {new MSPeak {XValue = shiftedMz, Height = 1}};
							for (double k = 1; k < 4; k++)
							{
								theoreticalPeakList.Add(new MSPeak {XValue = shiftedMz + (k*1.003/feature.ChargeState), Height = (float) (1.0 - (k/4))});
								theoreticalPeakList.Add(new MSPeak {XValue = shiftedMz - (k*1.003/feature.ChargeState), Height = (float) (1.0 - (k/4))});
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
							IsotopicProfile foundProfile = msFeatureFinder.FindMSFeature(candidatePeaks, isotopicProfile, 20, false);

							/*
							 * It is possible that the set mono pass of the previous theoretical distribution was the right-most peak of the actual distribution
							 * If so, we should be able to shift the theoretical distribution over to the left and find the actual distribution
							 */
							if (foundProfile == null)
							{
								foreach (MSPeak msPeak in sortPeaksQuery)
								{
									msPeak.XValue -= (1.003/feature.ChargeState);
								}

								isotopicProfile = new IsotopicProfile
								                  	{
								                  		MonoIsotopicMass = shiftedMass - 1.003,
								                  		MonoPeakMZ = shiftedMz - (1.003/feature.ChargeState),
								                  		ChargeState = feature.ChargeState,
								                  		Peaklist = sortPeaksQuery.ToList()
								                  	};

								foundProfile = msFeatureFinder.FindMSFeature(candidatePeaks, isotopicProfile, 20, false);
							}

							// Add to results, even if we did not find it.
							bool didFindProfile = foundProfile != null;
							crossLinkResult.MassShiftResults.KvpList.Add(new KeyValuePair<double, bool>(shiftedMass, didFindProfile));
						}

						crossLinkResultList.Add(crossLinkResult);
					}
				}
			}

			OutputCrossLinkResults(crossLinkResultList);
		}

		/// <summary>
		/// Use binary search to find all Isotopic Peaks we want to consider when looking for the cross-links shifts.
		/// </summary>
		/// <param name="completePeakList">The complete list of Isotopic Peaks that we will use for searching.</param>
		/// <param name="minimumMz">The minimum m/z value to consider.</param>
		/// <param name="scanLc">The LC Scan to consider.</param>
		/// <param name="scanIms">The IMS Scan to consider.</param>
		/// <returns>All peaks that are in the given LC Scan and IMS Scan and m/z >= thegiven m/z of the Feature.</returns>
		private static List<IPeak> FindCandidatePeaks(List<MSPeakResult> completePeakList, double minimumMz, int scanLc, int scanIms)
		{
			// Set up Peak Comparer to use for binary search later on
			AnonymousComparer<MSPeakResult> peakComparer = new AnonymousComparer<MSPeakResult>((x, y) => x.Frame_num != y.Frame_num ? x.Frame_num.CompareTo(y.Frame_num) : x.Scan_num != y.Scan_num ? x.Scan_num.CompareTo(y.Scan_num) : x.XValue.CompareTo(y.XValue));

			MSPeak msPeakLow = new MSPeak(minimumMz, 1, 1, 1);
			MSPeak msPeakHigh = new MSPeak(0, 1, 1, 1);
			MSPeakResult lowPeak = new MSPeakResult(1, scanLc, scanIms, msPeakLow);
			MSPeakResult highPeak = new MSPeakResult(1, scanLc, scanIms + 1, msPeakHigh);

			int lowPeakPosition = completePeakList.BinarySearch(lowPeak, peakComparer);
			int highPeakPosition = completePeakList.BinarySearch(highPeak, peakComparer);

			lowPeakPosition = lowPeakPosition < 0 ? ~lowPeakPosition : lowPeakPosition;
			highPeakPosition = highPeakPosition < 0 ? ~highPeakPosition : highPeakPosition;

			List<IPeak> candidatePeaks = new List<IPeak>();

			for (int j = lowPeakPosition; j < highPeakPosition; j++)
			{
				MSPeakResult msPeakResult = completePeakList[j];
				MSPeak msPeak = new MSPeak(msPeakResult.XValue, msPeakResult.Height, msPeakResult.Width, 1);
				candidatePeaks.Add(msPeak);
			}

			return candidatePeaks;
		} 

		/// <summary>
		/// Writes the results of cross-link searching to a csv file.
		/// </summary>
		/// <param name="crossLinkResultEnumerable">The List of CrossLinkResults objects</param>
		private static void OutputCrossLinkResults(IEnumerable<CrossLinkResult> crossLinkResultEnumerable)
		{
			TextWriter crossLinkWriter = new StreamWriter("crossLinkResults.csv");
			crossLinkWriter.WriteLine("Index,Pep1,Pep2,ModType,TheoreticalMass,FeatureMass,FeatureMz,ShiftedMassPep1,ShiftedMzPep1,ShiftedMassPep2,ShiftedMzPep2,ShiftedMassBoth,ShiftedMzBoth,ChargeState,LCScans,IMSScan,DriftTime,Abundance,FeatureIndex");
			int index = 0;

			var groupByCrossLinkQuery = crossLinkResultEnumerable.GroupBy(crossLinkResult => new { crossLinkResult.CrossLink, 
																								   crossLinkResult.LcImsMsFeature.FeatureId,
																								   crossLinkResult.MassShiftResults })
																 .Select(group => group.OrderBy(crossLinkResult => crossLinkResult.ScanLc) );

			foreach (var group in groupByCrossLinkQuery)
			{
				CrossLinkResult crossLinkResult = group.First();
				CrossLink crossLink = crossLinkResult.CrossLink;
				LcImsMsFeature feature = crossLinkResult.LcImsMsFeature;

				double featureMass = feature.MassMonoisotopic;
				double featureMz = (featureMass / feature.ChargeState) + MASS_OF_PROTON;

				StringBuilder outputLine = new StringBuilder();
				outputLine.Append(index++ + ",");
				outputLine.Append(crossLink.PeptideOne.SequenceOneLetter + ",");
				outputLine.Append((crossLink.PeptideTwo != null ? crossLink.PeptideTwo.SequenceOneLetter : "null") + ",");
				outputLine.Append(crossLink.ModType + ",");
				outputLine.Append(crossLink.Mass + ",");
				outputLine.Append(featureMass + ",");
				outputLine.Append(featureMz + ",");

				// Iterate over the results for each mass-shift
				foreach (KeyValuePair<double, bool> massShiftResult in crossLinkResult.MassShiftResults.KvpList)
				{
					// If the mass-shift value was found
					if (massShiftResult.Value)
					{
						double shiftMass = massShiftResult.Key;
						double shiftMz = (shiftMass / feature.ChargeState) + MASS_OF_PROTON;

						// Output information about the mass-shift
						outputLine.Append(shiftMass + "," + shiftMz + ",");
					}
					else
					{
						outputLine.Append("0,0,");
					}

					// If only 1 mass-shift, then the other 2 are N/A
					if (crossLinkResult.MassShiftResults.KvpList.Count == 1)
					{
						outputLine.Append("N/A,N/A,N/A,N/A,");
					}
				}

				outputLine.Append(feature.ChargeState + ",");

				foreach (CrossLinkResult individualCrossLinkResult in group)
				{
					outputLine.Append(individualCrossLinkResult.ScanLc + ";");
				}
				outputLine.Append(",");

				outputLine.Append(feature.ScanImsRep + ",");
				outputLine.Append(feature.DriftTime + ",");
				outputLine.Append(feature.Abundance + ",");
				outputLine.Append(feature.FeatureId + "\n");

				crossLinkWriter.Write(outputLine);
			}

			crossLinkWriter.Close();
		}
	}
}
