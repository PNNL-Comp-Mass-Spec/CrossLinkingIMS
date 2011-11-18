using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.IO;
using System.Linq;
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

			var sortPeakListQuery = from peak in peakList
									orderby peak.Frame_num, peak.Scan_num, peak.XValue
									select peak;

			peakList = sortPeakListQuery.ToList();

			// Set up Peak Comparer to use for binary search later on
			AnonymousComparer<MSPeakResult> peakComparer = new AnonymousComparer<MSPeakResult>((x, y) => x.Frame_num != y.Frame_num ? x.Frame_num.CompareTo(y.Frame_num) : x.Scan_num != y.Scan_num ? x.Scan_num.CompareTo(y.Scan_num) : x.XValue.CompareTo(y.XValue));

			TextWriter crossLinkWriter = new StreamWriter("crossLinkResults.csv");
			crossLinkWriter.WriteLine("Index,Pep1,Pep2,ModType,TheoreticalMass,FeatureMass,ShiftedMassPep1,ShiftedMzPep1,ShiftedMassPep2,ShiftedMzPep2,ShiftedMassBoth,ShiftedMzBoth,ChargeState,LCScan,IMSScan,DriftTime,Abundance,FeatureIndex");

			int index = 0;

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

					for (int currentScanLc = feature.ScanLcStart; currentScanLc <= feature.ScanLcEnd; currentScanLc++)
					{

						Console.WriteLine(crossLink.Mass + "\t" + feature.MassMonoisotopic + "\t" + feature.MzMonoisotopic + "\t" +
										  currentScanLc + "\t" + feature.ChargeState + "\t" + feature.DriftTime);

						List<IPeak> candidatePeaks = new List<IPeak>();

						/*
						 * Setup binary search to find all Isotopic Peaks we want to consider when looking for the cross-links shifts
						 *	- All peaks that are in the same LC Scan and IMS Scan as the Feature and m/z >= the monoisotopic m/x of the Feature
						 */
						MSPeak msPeakLow = new MSPeak(feature.MzMonoisotopic, 1, 1, 1);
						MSPeak msPeakHigh = new MSPeak(0, 1, 1, 1);
						MSPeakResult lowPeak = new MSPeakResult(1, currentScanLc, feature.ScanImsRep, msPeakLow);
						MSPeakResult highPeak = new MSPeakResult(1, currentScanLc, feature.ScanImsRep + 1, msPeakHigh);

						int lowPeakPosition = peakList.BinarySearch(lowPeak, peakComparer);
						int highPeakPosition = peakList.BinarySearch(highPeak, peakComparer);

						lowPeakPosition = lowPeakPosition < 0 ? ~lowPeakPosition : lowPeakPosition;
						highPeakPosition = highPeakPosition < 0 ? ~highPeakPosition : highPeakPosition;

						Console.WriteLine("Peaks " + lowPeakPosition + " to " + highPeakPosition);

						for (int j = lowPeakPosition; j < highPeakPosition; j++)
						{
							MSPeakResult msPeakResult = peakList[j];
							MSPeak msPeak = new MSPeak(msPeakResult.XValue, msPeakResult.Height, msPeakResult.Width, 1);
							candidatePeaks.Add(msPeak);
						}

						List<double> massShiftList = crossLink.MassShiftList;
						List<double> shiftedMassList = new List<double>();

						// Calculate the shifted mass values that we want to search for
						switch (massShiftList.Count)
						{
							case 1:
								{
									double firstNewMass = feature.MassMonoisotopic + massShiftList[0];
									//double secondNewMass = feature.MassMonoisotopic + massShiftList[0] + massShiftList[0];

									shiftedMassList.Add(firstNewMass);
									//shiftedMassList.Add(secondNewMass);
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

						List<Tuple<double, double, bool>> resultList = new List<Tuple<double, double, bool>>();

						// Search for shifted mass values in Isotopic Peaks
						foreach (double shiftedMass in shiftedMassList)
						{
							double shiftedMz = (shiftedMass/feature.ChargeState) + MASS_OF_PROTON;

							// Create theoretical Isotopic Peaks that will later form a theoretical Isotopic Profile
							List<MSPeak> theoreticalPeakList = new List<MSPeak> {new MSPeak {XValue = shiftedMz, Height = 1}};
							for (double k = 1; k < 4; k++)
							{
								theoreticalPeakList.Add(new MSPeak
								                        	{XValue = shiftedMz + (k*1.003/feature.ChargeState), Height = (float) (1.0 - (k/4))});
								theoreticalPeakList.Add(new MSPeak
								                        	{XValue = shiftedMz - (k*1.003/feature.ChargeState), Height = (float) (1.0 - (k/4))});
							}

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
							Console.WriteLine("Shifted Mass: " + shiftedMass + "\t" + shiftedMz + ": " +
							                  (foundProfile != null ? "true" : "false"));

							// The object will only be null if the Isotopic Profile was not found
							if (foundProfile != null)
							{
								//crossLinkWriter.WriteLine(crossLink.PeptideOne.SequenceOneLetter + "," + (crossLink.PeptideTwo != null ? crossLink.PeptideTwo.SequenceOneLetter : "null") + "," +
								//                          crossLink.ModType + "," + crossLink.Mass + "," + feature.MassMonoisotopic + "," + feature.MzMonoisotopic + "," + shiftedMass + "," + 
								//                          shiftedMz + "," + feature.ScanLcRep + "," + feature.ScanImsRep + "," + feature.DriftTime + "," + feature.ChargeState);
							}

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
								Console.WriteLine("Shifted Mass: " + shiftedMass + "\t" + shiftedMz + ": " +
								                  (foundProfile != null ? "true" : "false"));
								if (foundProfile != null)
								{
									//crossLinkWriter.WriteLine(crossLink.PeptideOne.SequenceOneLetter + "," + (crossLink.PeptideTwo != null ? crossLink.PeptideTwo.SequenceOneLetter : "null") + "," +
									//                          crossLink.ModType + "," + crossLink.Mass + "," + feature.MassMonoisotopic + "," + feature.MzMonoisotopic + "," + shiftedMass + "," +
									//                          shiftedMz + "," + feature.ScanLcRep + "," + feature.ScanImsRep + "," + feature.DriftTime + "," + feature.ChargeState);
								}
							}

							resultList.Add(new Tuple<double, double, bool>(shiftedMass, shiftedMz, foundProfile != null));
						}

						crossLinkWriter.Write(index++ + "," + crossLink.PeptideOne.SequenceOneLetter + "," +
						                      (crossLink.PeptideTwo != null ? crossLink.PeptideTwo.SequenceOneLetter : "null") + "," +
						                      crossLink.ModType + "," + crossLink.Mass + "," + feature.MassMonoisotopic + ",");

						switch (massShiftList.Count)
						{
							case 1:
								{
									if (resultList[0].Item3)
									{
										crossLinkWriter.Write(resultList[0].Item1 + "," + resultList[0].Item2 + ",");
									}
									else
									{
										crossLinkWriter.Write("0,0,");
									}

									crossLinkWriter.Write("N/A,N/A,N/A,N/A,");
								}
								break;
							case 2:
								{
									foreach (var tuple in resultList)
									{
										if (tuple.Item3)
										{
											crossLinkWriter.Write(tuple.Item1 + "," + tuple.Item2 + ",");
										}
										else
										{
											crossLinkWriter.Write("0,0,");
										}
									}
								}
								break;
						}

						crossLinkWriter.Write(feature.ChargeState + "," + currentScanLc + "," + feature.ScanImsRep + "," +
						                      feature.DriftTime + "," + feature.Abundance + "," + feature.FeatureId + "\n");
					}
				}
			}

			crossLinkWriter.Close();
		}
	}
}
