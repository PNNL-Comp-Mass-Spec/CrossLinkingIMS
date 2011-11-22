using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Threading;
using CrossLinkingIMS.Control;
using CrossLinkingIMS.Data;
using CrossLinkingIMS.Util;
using PNNLOmics.Utilities.ConsoleUtil;

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

		static void Main(string[] args)
		{
			// Handles a bug in .NET console applications
			IntPtr handle = Process.GetCurrentProcess().MainWindowHandle;
			SetConsoleMode(handle, ENABLE_EXTENDED_FLAGS);

			// Get the version number
			Assembly assembly = Assembly.GetExecutingAssembly();
			string assemblyVersion = assembly.GetName().Version.ToString();

			Console.WriteLine("CrossLinkingIMS Console Application Version " + assemblyVersion);

			CommandLineUtil commandLineUtil = new CommandLineUtil();

			// Make sure we have parameters. If not, then show the user the proper syntax
			bool doParametersExist = commandLineUtil.ParseCommandLine();
			if (!doParametersExist)
			{
				ShowSyntax();
				return;
			}

			// Get the Feature File Location
			string featureFileLocation = "";
			if (!commandLineUtil.RetrieveValueForParameter("f", out featureFileLocation))
			{
				Console.WriteLine("-f switch is missing");
				ShowSyntax();
				return;
			}

			FileInfo featureFile = new FileInfo(featureFileLocation);

			// Get the Peaks File Location
			string peaksFileLocation = "";
			if (!commandLineUtil.RetrieveValueForParameter("p", out peaksFileLocation))
			{
				Console.WriteLine("-p switch is missing");
				ShowSyntax();
				return;
			}

			FileInfo peaksFile = new FileInfo(peaksFileLocation);

			// Get the FastA File Location
			string fastAFileLocation = "";
			if (!commandLineUtil.RetrieveValueForParameter("fasta", out fastAFileLocation))
			{
				Console.WriteLine("-fasta switch is missing");
				ShowSyntax();
				return;
			}

			FileInfo fastAFile = new FileInfo(fastAFileLocation);

			// Hard-coded Protein Sequence
			//List<string> proteinList = new List<string>();
			//proteinList.Add("AEQVSKQEISHFKLVKVGTINVSQSGGQISSPSDLREKLSELADAKGGKYYHIIAAREHGPNFEAVAEVYNDATKLEHHHHHH");

			// Get the PPM Mass Tolerance
			string massToleranceString = "";
			if (!commandLineUtil.RetrieveValueForParameter("ppm", out massToleranceString))
			{
				Console.WriteLine("-ppm switch is missing");
				ShowSyntax();
				return;
			}

			double massTolerance = double.Parse(massToleranceString);

			// Run the cross-linking application
			Console.WriteLine("Executing...");
			IEnumerable<CrossLinkResult> crossLinkResults = CrossLinkingImsController.Execute(massTolerance, fastAFile, featureFile, peaksFile);

			// Hard coded output file
			FileInfo outputFileInfo = new FileInfo("crossLinkResults.csv");

			// Output the results
			Console.WriteLine("Outputting results to " + outputFileInfo.FullName);
			CrossLinkUtil.OutputCrossLinkResults(crossLinkResults, outputFileInfo);
		}

		private static void ShowSyntax()
		{
			Console.WriteLine();
			Console.WriteLine("Program syntax:");
			Console.WriteLine(Path.GetFileName(Assembly.GetExecutingAssembly().Location));
			Console.WriteLine("CrossLinkingIMSConsole.exe -f [Features File] -p [Peaks File] -fasta [FastA file] -ppm Value [optional arguments]");
			Console.WriteLine();
			Console.WriteLine("*********REQUIRED ARGUMENTS ***********.");
			Console.WriteLine();
			Console.WriteLine(" -f: Features File. LC-IMS-MS Feature Finder Output. See README.");
			Console.WriteLine(" -p: Peaks File. DeconTools Output. See README");
			Console.WriteLine(" -fasta: FastA File. Contains all protein sequences to search.");
			Console.WriteLine(" -ppm [value] : Mass tolerance in ppm");
			Console.WriteLine();
			Console.WriteLine("*********OPTIONAL ARGUMENTS ***********.");
			Console.WriteLine();
			Console.WriteLine(" -debug : Display detailed debug messages during iteration. (NOT YET IMPLEMENTED)");
			Console.WriteLine("");
			Thread.Sleep(2000);
		}
	}
}
