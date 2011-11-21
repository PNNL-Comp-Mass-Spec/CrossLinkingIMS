using System.Collections.Generic;
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
		public static List<IPeak> FindCandidatePeaks(List<MSPeakResult> completePeakList, double minimumMz, int scanLc, int scanIms)
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
	}
}
