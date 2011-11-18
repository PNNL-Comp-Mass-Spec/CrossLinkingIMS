namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Stores the number of carbon, hydrogen, nitrogen, oxygen, sulfur, and iron for a single sequence.
	/// </summary>
	public class AminoAcid
	{
		/// <summary>The number of carbon elements in the amino acid sequence</summary>
		public int NumCarbon { get; set; }

		/// <summary>The number of hydrogen elements in the amino acid sequence</summary>
		public int NumHydrogen { get; set; }

		/// <summary>The number of nitrogen elements in the amino acid sequence</summary>
		public int NumNitrogen { get; set; }

		/// <summary>The number of oxygen elements in the amino acid sequence</summary>
		public int NumOxygen { get; set; }

		/// <summary>The number of sulfur elements in the amino acid sequence</summary>
		public int NumSulfur { get; set; }

		/// <summary>The number of iron elements in the amino acid sequence</summary>
		public int NumIron { get; set; }
	}
}
