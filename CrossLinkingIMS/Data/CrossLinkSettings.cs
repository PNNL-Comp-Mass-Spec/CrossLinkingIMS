using System;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Class for storing the configurable parameters to be used by the cross-linking software.
	/// </summary>
	public class CrossLinkSettings
	{
		/// <summary>
		/// The mass tolerance in ppm
		/// </summary>
		public double MassTolerance { get; private set; }

		/// <summary>
		/// The trypsin digestion rule (fully tryptic, partially tryptic, no rules)
		/// </summary>
		public clsInSilicoDigest.CleavageRuleConstants TrypticType { get; private set; }

		/// <summary>
		/// The maximum number of missed cleavages to consider during digestion
		/// </summary>
		public int MaxMissedCleavages { get; private set; }

	    
        /// <summary>
        /// When True, crosslinked peptides have C13 labeling
        /// </summary>
        public bool UseC13 { get; private set; }
        
        /// <summary>
        /// When True, crosslinked peptides have N15 labeling
        /// </summary>
        public bool UseN15 { get; private set; }

	    /// <summary>
	    /// Default constructor which accepts all setting parameters.
	    /// </summary>
	    /// <param name="massTolerance">The mass tolerance in ppm</param>
	    /// <param name="maxMissedCleavages">The maximum number of missed cleavages to consider during digestion</param>
	    /// <param name="trypticString">The trypsin digestion rule (fully tryptic, partially tryptic, no rules)</param>
        /// <param name="useC13"> When True, crosslinked peptides have C13 labeling</param>
        /// <param name="useN15"> When True, crosslinked peptides have N15 labeling</param>
	    public CrossLinkSettings(double massTolerance, int maxMissedCleavages, string trypticString, bool useC13, bool useN15)
		{
			this.MassTolerance = massTolerance;
			this.MaxMissedCleavages = maxMissedCleavages;

	        this.UseC13 = useC13;
	        this.UseN15 = useN15;

			if(trypticString == null || trypticString.Equals("full"))
			{
				Console.WriteLine("Using Fully tryptic peptides only. " + maxMissedCleavages + " max missed cleavages.");
				this.TrypticType = clsInSilicoDigest.CleavageRuleConstants.ConventionalTrypsin;
			}
			else if (trypticString.Equals("partial"))
			{
				Console.WriteLine("Using Partially and Fully tryptic peptides. " + maxMissedCleavages + " max missed cleavages.");
				this.TrypticType = clsInSilicoDigest.CleavageRuleConstants.EricPartialTrypsin;
			}
			else if (trypticString.Equals("none"))
			{
				Console.WriteLine("Using no digestion rules. " + maxMissedCleavages + " max missed cleavages.");
				this.TrypticType = clsInSilicoDigest.CleavageRuleConstants.NoRule;
			}
			else
			{
				Console.WriteLine("Using Fully tryptic peptides only. " + maxMissedCleavages + " max missed cleavages.");
				this.TrypticType = clsInSilicoDigest.CleavageRuleConstants.ConventionalTrypsin;
			}
		}
	}
}
