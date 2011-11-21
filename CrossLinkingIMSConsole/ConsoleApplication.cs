using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.InteropServices;
using System.Text;
using CrossLinkingIMS.Control;
using CrossLinkingIMS.Data;
using CrossLinkingIMS.Util;

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

			Console.WriteLine("Using CrossLinkingIMS Version " + assemblyVersion);

			// Hard-coded mass tolerance
			const double massTolerance = 20;

			// Hard-coded Protein Sequence
			List<string> proteinList = new List<string>();
			proteinList.Add("AEQVSKQEISHFKLVKVGTINVSQSGGQISSPSDLREKLSELADAKGGKYYHIIAAREHGPNFEAVAEVYNDATKLEHHHHHH");

			// Hard-coded features file
			FileInfo featureFile = new FileInfo("SrfN_NC50_IMS_7Sep11_Roc_11-05-30_LCMSFeatures.txt");

			// Hard coded peaks file
			FileInfo peaksFile = new FileInfo("SrfN_NC50_IMS_7Sep11_Roc_11-05-30_peaks.txt");

			// Run the cross-linking application
			Console.WriteLine("Executing...");
			IEnumerable<CrossLinkResult> crossLinkResults = CrossLinkingImsController.Execute(massTolerance, proteinList, featureFile, peaksFile);

			// Hard coded output file
			FileInfo outputFileInfo = new FileInfo("crossLinkResults.csv");

			// Output the results
			Console.WriteLine("Outputting results to " + outputFileInfo.FullName);
			CrossLinkUtil.OutputCrossLinkResults(crossLinkResults, outputFileInfo);
		}
	}
}
