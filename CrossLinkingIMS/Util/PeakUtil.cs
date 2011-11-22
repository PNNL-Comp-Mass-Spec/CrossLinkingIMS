using System.Collections.Generic;
using CrossLinkingIMS.Data;
using DeconTools.Backend.Core;
using DeconTools.Backend.DTO;

namespace CrossLinkingIMS.Util
{
	/// <summary>
	/// Utility class containing function related to Isotopic Peaks.
	/// </summary>
	public class PeakUtil
	{
		/// <summary>
		/// Use binary search to find all Isotopic Peaks we want to consider when looking for the cross-links shifts.
		/// </summary>
		/// <param name="completePeakList">The complete list of Isotopic Peaks that we will use for searching.</param>
		/// <param name="minimumMz">The minimum m/z value to consider.</param>
		/// <param name="scanLc">The LC Scan to consider.</param>
		/// <param name="scanIms">The IMS Scan to consider.</param>
		/// <returns>All peaks that are in the given LC Scan and IMS Scan and m/z >= thegiven m/z of the Feature.</returns>
		public static List<IPeak> FindCandidatePeaks(List<IsotopicPeak> completePeakList, double minimumMz, int scanLc, int scanIms)
		{
			// Set up Peak Comparer to use for binary search later on
			AnonymousComparer<IsotopicPeak> peakComparer = new AnonymousComparer<IsotopicPeak>((x, y) => x.ScanLc != y.ScanLc ? x.ScanLc.CompareTo(y.ScanLc) : x.ScanIms != y.ScanIms ? x.ScanIms.CompareTo(y.ScanIms) : x.Mz.CompareTo(y.Mz));

			IsotopicPeak lowPeak = new IsotopicPeak {ScanLc = scanLc, ScanIms = scanIms, Mz = minimumMz, Intensity = 1};
			IsotopicPeak highPeak = new IsotopicPeak { ScanLc = scanLc, ScanIms = scanIms + 1, Mz = 0, Intensity = 1 };

			int lowPeakPosition = completePeakList.BinarySearch(lowPeak, peakComparer);
			int highPeakPosition = completePeakList.BinarySearch(highPeak, peakComparer);

			lowPeakPosition = lowPeakPosition < 0 ? ~lowPeakPosition : lowPeakPosition;
			highPeakPosition = highPeakPosition < 0 ? ~highPeakPosition : highPeakPosition;

			List<IPeak> candidatePeaks = new List<IPeak>();

			for (int j = lowPeakPosition; j < highPeakPosition; j++)
			{
				IsotopicPeak peak = completePeakList[j];
				MSPeak msPeak = new MSPeak(peak.Mz, peak.Intensity, 0.05f, 1);
				candidatePeaks.Add(msPeak);
			}

			return candidatePeaks;
		}
	}
}
