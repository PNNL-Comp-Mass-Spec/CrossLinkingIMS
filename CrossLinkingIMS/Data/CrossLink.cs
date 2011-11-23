using System.Collections.Generic;
using CrossLinkingIMS.Util;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Specifies the type of Modification of the cross-link.
	/// See Figure 2 of http://www.springerlink.com/content/g31r736110816733/fulltext.pdf for more information.
	/// </summary>
	public enum ModType
	{
		/// <summary>Unmodified peptide</summary>
		None,

		/// <summary>One end of cross-linking reagent reacted with the protein, the other reacted with water (often called a “dead-end”).</summary>
		Zero,

		/// <summary>Each ends of the reagent reacted with the same peptide (intrapeptide cross-link).</summary>
		One,

		/// <summary>Each end of the reagent reacted with a different peptide (interpeptide cross-link).</summary>
		Two,

		/// <summary>A mix of Type 0 and Type 1 Modifications</summary>
		ZeroOne,

		/// <summary>A mix of Type 0 and Type 2 Modifications</summary>
		ZeroTwo
	};

	/// <summary>
	/// Stores information about a single cross-link for either 1 or 2 peptides.
	/// </summary>
	public class CrossLink
	{
		/// <summary>
		/// The identifier of the protein used to generate this cross-link
		/// </summary>
		public string ProteinId { get; private set; }

		/// <summary>
		/// The first Peptide of the cross-link.
		/// </summary>
		public clsInSilicoDigest.PeptideInfoClass PeptideOne { get; private set; }

		/// <summary>
		/// The second peptide of the cross-link. May be null.
		/// </summary>
		public clsInSilicoDigest.PeptideInfoClass PeptideTwo { get; private set; }

		/// <summary>
		/// The monoisotopic mass of the un-shifted cross-link, in daltons.
		/// </summary>
		public double Mass { get; private set; }

		/// <summary>
		/// The mod type of the cross link. See CrossLink.ModType for explaination of mod types.
		/// </summary>
		public ModType ModType { get; private set; }

		/// <summary>
		/// The List of possible shifts in mass for the cross-link to exhibit.
		/// </summary>
		public List<double> MassShiftList { get; private set; }

		/// <summary>
		/// The only constructor for the CrossLink object requires 1 or 2 peptides, the mass of the peptides, and the ModType of the cross-link.
		/// </summary>
		/// <param name="proteinId">The identifier of the protein of the cross-link.</param>
		/// <param name="peptideOne">The first peptide of the cross-link.</param>
		/// <param name="peptideTwo">The second peptide of the cross-link. Can be null if linking the first peptide to itself.</param>
		/// <param name="mass">The monoisotopic mass of the un-shifted cross-link, in daltons.</param>
		/// <param name="modType">The mod type of the cross link. See CrossLink.ModType for explaination of mod types.</param>
		public CrossLink(string proteinId, clsInSilicoDigest.PeptideInfoClass peptideOne, clsInSilicoDigest.PeptideInfoClass peptideTwo, double mass, ModType modType)
		{
			this.ProteinId = proteinId;
			this.PeptideOne = peptideOne;
			this.PeptideTwo = peptideTwo;
			this.Mass = mass;
			this.ModType = modType;

			this.MassShiftList = new List<double>();
			if (peptideOne != null) this.MassShiftList.Add(CrossLinkUtil.CalculateMassShift(peptideOne.SequenceOneLetter));
			if (peptideTwo != null) this.MassShiftList.Add(CrossLinkUtil.CalculateMassShift(peptideTwo.SequenceOneLetter));
		}

		/// <summary>
		/// Determines whether the specified Object is equal to the current Object.
		/// Equality is determined by Mass and ModType.
		/// </summary>
		/// <param name="obj">The specified Object.</param>
		/// <returns>true if the objects are equal, false otherwise</returns>
		public override bool Equals(object obj)
		{
			if (ReferenceEquals(null, obj)) return false;
			if (ReferenceEquals(this, obj)) return true;
			if (obj.GetType() != typeof (CrossLink)) return false;
			return Equals((CrossLink) obj);
		}

		/// <summary>
		/// Determines whether the specified CrossLink Object is equal to the current CrossLink Object.
		/// Equality is determined by Protein, Mass, and ModType.
		/// </summary>
		/// <param name="other">The specified CrossLink Object.</param>
		/// <returns>true if the objects are equal, false otherwise</returns>
		public bool Equals(CrossLink other)
		{
			if (ReferenceEquals(null, other)) return false;
			if (ReferenceEquals(this, other)) return true;
			return other.ProteinId.Equals(this.ProteinId) && other.Mass.Equals(this.Mass) && Equals(other.ModType, this.ModType);
		}

		public override int GetHashCode()
		{
			unchecked
			{
				return this.ProteinId.GetHashCode()*397 ^ this.Mass.GetHashCode()*397 ^ this.ModType.GetHashCode();
			}
		}

		public override string ToString()
		{
			return string.Format("ProteinId: {0}, PeptideOne: {1}, PeptideTwo: {2}, ModType: {3}", ProteinId, PeptideOne.SequenceOneLetter, PeptideTwo != null ? PeptideTwo.SequenceOneLetter : "", ModType);
		}
	}
}
