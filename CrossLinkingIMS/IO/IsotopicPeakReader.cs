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
        private const string IMS_SCAN_NUM = "ims_scan_num";
        private const string MZ = "mz";
        private const string INTENSITY = "intensity";
        private const string ABUNDANCE = "abundance";

        /// <summary>
        /// Reads a peaks file and returns the list of Isotopic Peaks.
        /// </summary>
        /// <param name="peakFile">The FileInfo object for the file to be read.</param>
        /// <returns>A List of IsotopicPeak objects read in from the file.</returns>
        public static List<IsotopicPeak> ReadFile(FileInfo peakFile)
        {
            var peakList = new List<IsotopicPeak>();

            using (TextReader textReader = new StreamReader(peakFile.FullName))
            {
                var columnHeaders = textReader.ReadLine();
                var columnMapping = CreateColumnMapping(columnHeaders);

                string line;
                while ((line = textReader.ReadLine()) != null)
                {
                    var peak = ParseLine(line, columnMapping);
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
        private static Dictionary<string, int> CreateColumnMapping(string columnString)
        {
            var columnMap = new Dictionary<string, int>();
            var columnTitles = columnString.Split('\t', ',', '\n');

            for (var i = 0; i < columnTitles.Count(); i++)
            {
                var columnTitle = columnTitles[i];

                switch (columnTitle)
                {
                    case FRAME_NUM:
                        columnMap.Add(FRAME_NUM, i);
                        break;
                    case SCAN_NUM:
                    case IMS_SCAN_NUM:
                        columnMap.Add(SCAN_NUM, i);
                        break;
                    case MZ:
                        columnMap.Add(MZ, i);
                        break;
                    case ABUNDANCE:
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
        private static IsotopicPeak ParseLine(string line, IDictionary<string, int> columnMapping)
        {
            var columnName = "unknown";
                
            try
            {
                var columns = line.Split('\t', ',', '\n');

                var peak = new IsotopicPeak();

                columnName = "FrameNum";
                if (columnMapping.ContainsKey(FRAME_NUM)) 
                    peak.ScanLc = int.Parse(columns[columnMapping[FRAME_NUM]]);

                columnName = "ScanNum";
                if (columnMapping.ContainsKey(SCAN_NUM)) 
                    peak.ScanIms = int.Parse(columns[columnMapping[SCAN_NUM]]);

                columnName = "MZ";
                if (columnMapping.ContainsKey(MZ)) 
                    peak.Mz = double.Parse(columns[columnMapping[MZ]]);

                columnName = "Intensity";
                if (columnMapping.ContainsKey(INTENSITY)) 
                    peak.Intensity = (int)Math.Round(double.Parse(columns[columnMapping[INTENSITY]]));

                return peak;

            }
            catch (Exception)
            {
                Console.WriteLine("Error parsing column " + columnName + " in line: " + line);
                throw;
            }

        }
    }
}
