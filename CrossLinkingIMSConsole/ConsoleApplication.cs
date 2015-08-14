using System;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Threading;
using CrossLinkingIMS.Control;
using CrossLinkingIMS.Data;
using CrossLinkingIMS.Util;
using PNNLOmics.Utilities;

namespace CrossLinkingIMSConsole
{
    /// <summary>
    /// Console Application for the CrossLinkingIMS DLL.
    /// </summary>
    class ConsoleApplication
    {
        [DllImport("kernel32.dll")]
        public static extern bool SetConsoleMode(IntPtr hConsoleHandle, uint dwMode);

        private const uint ENABLE_EXTENDED_FLAGS = 0x0080;

        private const string PROGRAM_DATE = "August 14, 2015";

        static int Main(string[] args)
        {
            try
            {

                // Handles a bug in .NET console applications
                var handle = Process.GetCurrentProcess().MainWindowHandle;
                SetConsoleMode(handle, ENABLE_EXTENDED_FLAGS);

                // Get the version number
                var assembly = Assembly.GetExecutingAssembly();
                var assemblyVersion = assembly.GetName().Version.ToString();

                Console.WriteLine("CrossLinkingIMS Console Application Version " + assemblyVersion);

                var commandLineUtil = new CommandLineUtil(caseSensitiveSwitchNames: false);

                // Make sure we have parameters. If not, then show the user the proper syntax
                var doParametersExist = commandLineUtil.ParseCommandLine();
                if (!doParametersExist)
                {
                    ShowSyntax();
                    return -1;
                }

                // Get the Feature File Location
                string featureFileLocation;
                if (!commandLineUtil.RetrieveValueForParameter("f", out featureFileLocation))
                {
                    Console.WriteLine("-f switch is missing");
                    ShowSyntax();
                    return -2;
                }

                if (!FileExists(featureFileLocation))
                {
                    return -4;
                }

                var featureFile = new FileInfo(featureFileLocation);

                // Get the Peaks File Location
                string peaksFileLocation;
                if (!commandLineUtil.RetrieveValueForParameter("p", out peaksFileLocation))
                {
                    Console.WriteLine("-p switch is missing");
                    ShowSyntax();
                    return -2;
                }

                if (!FileExists(peaksFileLocation))
                {
                    return -4;
                }

                var peaksFile = new FileInfo(peaksFileLocation);

                // Get the FastA File Location
                string fastAFileLocation;
                if (!commandLineUtil.RetrieveValueForParameter("fasta", out fastAFileLocation))
                {
                    Console.WriteLine("-fasta switch is missing");
                    ShowSyntax();
                    return -2;
                }

                if (!FileExists(fastAFileLocation))
                {
                    return -4;
                }

                var fastAFile = new FileInfo(fastAFileLocation);

                // Get the PPM Mass Tolerance
                string massToleranceString;
                if (!commandLineUtil.RetrieveValueForParameter("ppm", out massToleranceString))
                {
                    Console.WriteLine("-ppm switch is missing");
                    ShowSyntax();
                    return -2;
                }

                var massTolerance = double.Parse(massToleranceString);

                // Get the Max Missed Cleavages
                string maxMissedCleavagesString;
                commandLineUtil.RetrieveValueForParameter("c", out maxMissedCleavagesString);

                if (string.IsNullOrEmpty(maxMissedCleavagesString))
                    maxMissedCleavagesString = "1";

                var maxMissedCleavages = int.Parse(maxMissedCleavagesString);

                // Get the Partially Tryptic Flag
                string trypticString;
                commandLineUtil.RetrieveValueForParameter("t", out trypticString);

                string staticDeltaMassString;
                var staticDeltaMass = 0.0;
                if (commandLineUtil.RetrieveValueForParameter("static", out staticDeltaMassString))
                {
                    double.TryParse(staticDeltaMassString, out staticDeltaMass);
                }

                var useC13 = true;
                var useN15 = true;
            
                if (commandLineUtil.IsParameterPresent("c13off"))
                    useC13 = false;

                if (commandLineUtil.IsParameterPresent("n15off"))
                    useN15 = false;

                if (!(useC13 || useN15) && Math.Abs(staticDeltaMass) < float.Epsilon)
                {
                    Console.WriteLine("If you use both -C13off and -N15off, you must use -Static:MassDiff to specify a static mass difference");
                    Thread.Sleep(1500);
                    return -3;
                }

                var settings = new CrossLinkSettings(massTolerance, maxMissedCleavages, trypticString, useC13, useN15, staticDeltaMass);

                // Get the Output File Location
                string outputFileLocation;
                if (!commandLineUtil.RetrieveValueForParameter("o", out outputFileLocation))
                {
                    if (featureFile.DirectoryName != null && featureFile.DirectoryName == peaksFile.DirectoryName)
                        outputFileLocation = Path.Combine(featureFile.DirectoryName, "crossLinkResults.csv");
                    else
                        outputFileLocation = "crossLinkResults.csv";
                }

                // Run the cross-linking application
                Console.WriteLine("Executing...");
                var crossLinkResults = CrossLinkingImsController.Execute(settings, fastAFile, featureFile, peaksFile);

                var outputFileInfo = new FileInfo(outputFileLocation);

                // Output the results
                Console.WriteLine("Outputting " + crossLinkResults.Count().ToString("#,##0") + " results to \n" + outputFileInfo.FullName);
                CrossLinkUtil.OutputCrossLinkResults(crossLinkResults, outputFileInfo);

                return 0;
            }
            catch (Exception ex)
            {
                Console.WriteLine();
                Console.WriteLine("=========================================");
                Console.WriteLine("Error: " + ex.Message);
                Console.WriteLine(ex.StackTrace);
                Console.WriteLine("=========================================");

                Thread.Sleep(2000);
                return -10;
            }

        }

        private static bool FileExists(string filePath)
        {
            if (File.Exists(filePath))
            {
                return true;
            }

            Console.WriteLine("File not found: " + filePath);
            Thread.Sleep(1500);

            return false;
        }

        private static string GetAppVersion()
        {
            return Assembly.GetExecutingAssembly().GetName().Version + " (" + PROGRAM_DATE + ")";
        }

        private static void ShowSyntax()
        {
            var exeName = Path.GetFileName(Assembly.GetExecutingAssembly().Location);

            Console.WriteLine();
            Console.WriteLine("Program syntax:");
            Console.WriteLine(exeName + "-f [Features File] -p [Peaks File] -fasta [FastA file] -ppm Value [optional arguments]");
            Console.WriteLine();
            Console.WriteLine("*********REQUIRED ARGUMENTS ***********");
            Console.WriteLine();
            Console.WriteLine(" -f        Features File. LC-IMS-MS Feature Finder Output. See README.");
            Console.WriteLine(" -p        Peaks File. DeconTools Output. See README");
            Console.WriteLine(" -fasta    FastA File. Contains all protein sequences to search.");
            Console.WriteLine(" -ppm [value]   Mass tolerance in ppm");
            Console.WriteLine();
            Console.WriteLine("*********OPTIONAL ARGUMENTS ***********");
            Console.WriteLine();
            Console.WriteLine(" -c   The maximum number of missed cleavages to consider. Defaults to 1.");
            Console.WriteLine(" -t   Set to 'full' for fully tryptic only. ");
            Console.WriteLine("      Set to 'partial' to consider partially and fully tryptic.");
            Console.WriteLine("      Set to 'none' to consider non, partially, and fully tryptic.");
            Console.WriteLine("      Defaults to 'full'.");
            Console.WriteLine(" -o   The desired location and name for the output file. Defaults to workingDirectory/crossLinkResults.csv");
            Console.WriteLine();
            Console.WriteLine("Crosslinked peptide pairs are found using a delta mass value that is");
            Console.WriteLine("typically computed based on the number of carbons and number of nitrogen atoms,");
            Console.WriteLine("treating carbons in the crosslinked peptide as C13 and nitrogens as N15.");
            Console.WriteLine("Alter this behavior using switches -C13off, -N15off, and -Static"); 
            Console.WriteLine(" -C13off           Peptides are not labeled with C13");
            Console.WriteLine(" -N15off           Peptides are not labeled with N15");
            Console.WriteLine(" -Static [Value]   Static (fixed) mass shift to use for all paired peptides");
            Console.WriteLine();
            Console.WriteLine("For example, to disable C13 and N15 shifts and only use a ");
            Console.WriteLine("fixed shift of 25.5 Da, your command line would look like this:");
            Console.WriteLine();
            Console.WriteLine(exeName);
            Console.WriteLine(" -f Dataset_LCMSFeatures.txt -p Dataset_Filtered_isos.csv -fasta Database.fasta");
            Console.WriteLine(" -ppm 10 -C13Off -N15off -Static 25.5");
            Console.WriteLine();
            Console.WriteLine("Program written by Kevin Crowell for the Department of Energy (PNNL, Richland, WA) in 2011");
            Console.WriteLine("Version: " + GetAppVersion());

            Console.WriteLine();

            Console.WriteLine("E-mail: matthew.monroe@pnnl.gov or matt@alchemistmatt.com");
            Console.WriteLine("Website: http://omics.pnl.gov/ or http://panomics.pnnl.gov/");
            Console.WriteLine();

            Thread.Sleep(1500);
        }
    }
}
