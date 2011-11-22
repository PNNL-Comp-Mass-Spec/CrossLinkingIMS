using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Stores information about a single Isotopic Peak
	/// </summary>
	public class IsotopicPeak
	{
		public int ScanLc { get; set; }
		public int ScanIms { get; set; }
		public double Mz { get; set; }
		public int Intensity { get; set; }
	}
}
