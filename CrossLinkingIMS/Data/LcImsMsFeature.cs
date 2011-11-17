using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Stores information about a single LC-IMS-MS Feature.
	/// </summary>
	public class LcImsMsFeature
	{
		public int FeatureId { get; set; }
		public int ChargeState { get; set; }
		public int ScanLcStart { get; set; }
		public int ScanLcEnd { get; set; }
		public int ScanLcRep { get; set; }
		public int ScanImsRep { get; set; }
		public double MassMonoisotopic { get; set; }
		public double MzMonoisotopic { get; set; }
		public double DriftTime { get; set; }
		public double Abundance { get; set; }
	}
}
