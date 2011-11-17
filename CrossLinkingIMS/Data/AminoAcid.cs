namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Stores the number of carbon, hydrogen, nitrogen, oxygen, sulfur, and iron for a single sequence.
	/// </summary>
	public class AminoAcid
	{
		public int NumCarbon { get; set; }
		public int NumHydrogen { get; set; }
		public int NumNitrogen { get; set; }
		public int NumOxygen { get; set; }
		public int NumSulfur { get; set; }
		public int NumIron { get; set; }
	}
}
