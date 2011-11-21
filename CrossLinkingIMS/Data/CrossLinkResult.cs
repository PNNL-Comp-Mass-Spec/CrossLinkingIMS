namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Stores information abotu the results of searching for a CrossLink in the data.
	/// </summary>
	public class CrossLinkResult
	{
		/// <summary>
		/// The CrossLink object that is being searched for.
		/// </summary>
		public CrossLink CrossLink { get; private set; }

		/// <summary>
		/// The LC-IMS-MS Feature that was found to be the intiial marker of the CrossLink.
		/// </summary>
		public LcImsMsFeature LcImsMsFeature { get; private set; }

		/// <summary>
		/// The results showing the mass shifts that were searched for and if they were found.
		/// </summary>
		public MassShiftResult MassShiftResults { get; set; }

		/// <summary>
		/// The specific LC Scan for these results.
		/// </summary>
		public int ScanLc { get; private set; }

		/// <summary>
		/// This only constructor for CrossLinkResult requires the CrossLink, LC-IMS-MS Feature, and LC Scan.
		/// The Mass Shift Results are added later on.
		/// </summary>
		/// <param name="crossLink"></param>
		/// <param name="lcImsMsFeature"></param>
		/// <param name="scanLc"></param>
		public CrossLinkResult(CrossLink crossLink, LcImsMsFeature lcImsMsFeature, int scanLc)
		{
			this.CrossLink = crossLink;
			this.LcImsMsFeature = lcImsMsFeature;
			this.ScanLc = scanLc;

			this.MassShiftResults = new MassShiftResult();
		}
	}
}
