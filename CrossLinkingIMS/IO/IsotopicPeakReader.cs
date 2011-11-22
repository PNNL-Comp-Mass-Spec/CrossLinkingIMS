using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using CrossLinkingIMS.Data;

namespace CrossLinkingIMS.IO
{
	/// <summary>
	/// Used for reading a text file that contains information about isotopic Peaks and creating IsotopicPeak objects.
	/// </summary>
	public class IsotopicPeakReader
	{
		private const string FRAME_NUM = "frame_num";
		private const string SCAN_NUM = "scan_num";
		private const string MZ = "mz";
		private const string INTENSITY = "intensity";

		/// <summary>
		/// Reads a peaks file and returns the list of Isotopic Peaks.
		/// </summary>
		/// <param name="peakFile">The FileInfo object for the file to be read.</param>
		/// <returns>A List of IsotopicPeak objects read in from the file.</returns>
		public static List<IsotopicPeak> ReadFile(FileInfo peakFile)
		{
			List<IsotopicPeak> peakList = new List<IsotopicPeak>();

			using (TextReader textReader = new StreamReader(peakFile.FullName))
			{
				string columnHeaders = textReader.ReadLine();
				Dictionary<string, int> columnMapping = CreateColumnMapping(columnHeaders);

				string line = "";
				while ((line = textReader.ReadLine()) != null)
				{
					IsotopicPeak peak = ParseLine(line, columnMapping);
					peakList.Add(peak);
				}
			}

			return peakList;
		}

		/// <summary>
		/// Using the header row of the input file, creates a mapping of column headers to their position in the file.
		/// </summary>
		/// <param name="columnString">The header row of the input file, as a string.</param>
		/// <returns>A Dictionary object containing the mapping between the headers and their column number position.</returns>
		private static Dictionary<string, int> CreateColumnMapping(String columnString)
		{
			Dictionary<string, int> columnMap = new Dictionary<string, int>();
			string[] columnTitles = columnString.Split('\t', ',', '\n');

			for (int i = 0; i < columnTitles.Count(); i++)
			{
				String columnTitle = columnTitles[i];

				switch (columnTitle)
				{
					case FRAME_NUM:
						columnMap.Add(FRAME_NUM, i);
						break;
					case SCAN_NUM:
						columnMap.Add(SCAN_NUM, i);
						break;
					case MZ:
						columnMap.Add(MZ, i);
						break;
					case INTENSITY:
						columnMap.Add(INTENSITY, i);
						break;
				}
			}

			return columnMap;
		}

		/// <summary>
		/// Parses a single line of the input file and returns a single Isotopic Peak.
		/// </summary>
		/// <param name="line">The line to read, as a string.</param>
		/// <param name="columnMapping">The column mapping for the input file.</param>
		/// <returns>A single IsotopicPeak object created by using the information from the parsed line.</returns>
		private static IsotopicPeak ParseLine(String line, IDictionary<string, int> columnMapping)
		{
			string[] columns = line.Split('\t');

			IsotopicPeak peak = new IsotopicPeak();

			if (columnMapping.ContainsKey(FRAME_NUM)) peak.ScanLc = int.Parse(columns[columnMapping[FRAME_NUM]]);
			if (columnMapping.ContainsKey(SCAN_NUM)) peak.ScanIms = int.Parse(columns[columnMapping[SCAN_NUM]]);
			if (columnMapping.ContainsKey(MZ)) peak.Mz = double.Parse(columns[columnMapping[MZ]]);
			if (columnMapping.ContainsKey(INTENSITY)) peak.Intensity = int.Parse(columns[columnMapping[INTENSITY]]);

			return peak;
		}
	}
}
