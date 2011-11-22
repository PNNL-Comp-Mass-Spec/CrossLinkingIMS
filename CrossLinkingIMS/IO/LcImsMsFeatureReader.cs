using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using CrossLinkingIMS.Constants;
using CrossLinkingIMS.Data;

namespace CrossLinkingIMS.IO
{
	/// <summary>
	/// Used for reading a text file that contains information about LC-IMS-MS Features and creating LC-IMS-MS Feature objects.
	/// </summary>
	public class LcImsMsFeatureReader
	{
		private const string FEATURE_INDEX = "Feature_Index";
		private const string MONOISOTOPIC_MASS = "Monoisotopic_Mass";
		private const string LC_SCAN_START = "Scan_Start";
		private const string LC_SCAN_END = "Scan_End";
		private const string LC_SCAN_REP = "Scan";
		private const string IMS_SCAN_REP = "IMS_Scan";
		private const string CHARGE_STATE = "Class_Rep_Charge";
		private const string DRIFT_TIME = "Drift_Time";
		private const string ABUNDANCE = "Abundance";

		/// <summary>
		/// Reads an LC-IMS-MS Feature file and returns the list of LC-IMS-MS Features.
		/// </summary>
		/// <param name="featureFile">The FileInfo object for the file to be read.</param>
		/// <returns>A List of LcImsMsFeature objects read in from the file.</returns>
		public static List<LcImsMsFeature> ReadFile(FileInfo featureFile)
		{
			List<LcImsMsFeature> featureList = new List<LcImsMsFeature>();

			using (TextReader textReader = new StreamReader(featureFile.FullName))
			{
				string columnHeaders = textReader.ReadLine();
				Dictionary<string, int> columnMapping = CreateColumnMapping(columnHeaders);

				string line = "";
				while ((line = textReader.ReadLine()) != null)
				{
					LcImsMsFeature feature = ParseLine(line, columnMapping);
					featureList.Add(feature);
				}
			}

			return featureList;
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
					case FEATURE_INDEX:
						columnMap.Add(FEATURE_INDEX, i);
						break;
					case MONOISOTOPIC_MASS:
						columnMap.Add(MONOISOTOPIC_MASS, i);
						break;
					case LC_SCAN_START:
						columnMap.Add(LC_SCAN_START, i);
						break;
					case LC_SCAN_END:
						columnMap.Add(LC_SCAN_END, i);
						break;
					case LC_SCAN_REP:
						columnMap.Add(LC_SCAN_REP, i);
						break;
					case IMS_SCAN_REP:
						columnMap.Add(IMS_SCAN_REP, i);
						break;
					case CHARGE_STATE:
						columnMap.Add(CHARGE_STATE, i);
						break;
					case DRIFT_TIME:
						columnMap.Add(DRIFT_TIME, i);
						break;
					case ABUNDANCE:
						columnMap.Add(ABUNDANCE, i);
						break;
				}
			}

			return columnMap;
		}

		/// <summary>
		/// Parses a single line of the input file and returns a single LC-IMS-MS Feature.
		/// </summary>
		/// <param name="line">The line to read, as a string.</param>
		/// <param name="columnMapping">The column mapping for the input file.</param>
		/// <returns>A single LcImsMsFeature object created by using the information from the parsed line.</returns>
		private static LcImsMsFeature ParseLine(String line, IDictionary<string, int> columnMapping)
		{
			string[] columns = line.Split('\t');

			LcImsMsFeature lcImsMsFeature = new LcImsMsFeature();

			if (columnMapping.ContainsKey(FEATURE_INDEX)) lcImsMsFeature.FeatureId = int.Parse(columns[columnMapping[FEATURE_INDEX]]);
			if (columnMapping.ContainsKey(MONOISOTOPIC_MASS)) lcImsMsFeature.MassMonoisotopic = double.Parse(columns[columnMapping[MONOISOTOPIC_MASS]]);
			if (columnMapping.ContainsKey(LC_SCAN_START)) lcImsMsFeature.ScanLcStart = int.Parse(columns[columnMapping[LC_SCAN_START]]);
			if (columnMapping.ContainsKey(LC_SCAN_END)) lcImsMsFeature.ScanLcEnd = int.Parse(columns[columnMapping[LC_SCAN_END]]);
			if (columnMapping.ContainsKey(LC_SCAN_REP)) lcImsMsFeature.ScanLcRep = int.Parse(columns[columnMapping[LC_SCAN_REP]]);
			if (columnMapping.ContainsKey(IMS_SCAN_REP)) lcImsMsFeature.ScanImsRep = int.Parse(columns[columnMapping[IMS_SCAN_REP]]);
			if (columnMapping.ContainsKey(CHARGE_STATE)) lcImsMsFeature.ChargeState = int.Parse(columns[columnMapping[CHARGE_STATE]]);
			if (columnMapping.ContainsKey(DRIFT_TIME)) lcImsMsFeature.DriftTime = double.Parse(columns[columnMapping[DRIFT_TIME]]);
			if (columnMapping.ContainsKey(ABUNDANCE)) lcImsMsFeature.Abundance = double.Parse(columns[columnMapping[ABUNDANCE]]);

			if (columnMapping.ContainsKey(MONOISOTOPIC_MASS) && columnMapping.ContainsKey(CHARGE_STATE))
			{
				double mz = (lcImsMsFeature.MassMonoisotopic / lcImsMsFeature.ChargeState) + GeneralConstants.MASS_OF_PROTON;
				lcImsMsFeature.MzMonoisotopic = mz;
			}

			return lcImsMsFeature;
		}
	}
}
