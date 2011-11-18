using System.Collections.Generic;
using CrossLinkingIMS.Util;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Specifies the type of Modification of the cross-link
	/// </summary>
	public enum ModType
	{
		/// <summary>Unmodified peptide</summary>
		None,

		/// <summary></summary>
		Zero, 

		/// <summary></summary>
		One,

		/// <summary></summary>
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
		/// <param name="peptideOne">The first peptide of the cross-link.</param>
		/// <param name="peptideTwo">The second peptide of the cross-link. Can be null if linking the first peptide to itself.</param>
		/// <param name="mass">The monoisotopic mass of the un-shifted cross-link, in daltons.</param>
		/// <param name="modType">The mod type of the cross link. See CrossLink.ModType for explaination of mod types.</param>
		public CrossLink(clsInSilicoDigest.PeptideInfoClass peptideOne, clsInSilicoDigest.PeptideInfoClass peptideTwo, double mass, ModType modType)
		{
			this.PeptideOne = peptideOne;
			this.PeptideTwo = peptideTwo;
			this.Mass = mass;
			this.ModType = modType;

			this.MassShiftList = new List<double>();
			if (peptideOne != null) this.MassShiftList.Add(SequenceUtil.CalculateMassShift(peptideOne.SequenceOneLetter));
			if (peptideTwo != null) this.MassShiftList.Add(SequenceUtil.CalculateMassShift(peptideTwo.SequenceOneLetter));
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
		/// Equality is determined by Mass and ModType.
		/// </summary>
		/// <param name="other">The specified CrossLink Object.</param>
		/// <returns>true if the objects are equal, false otherwise</returns>
		public bool Equals(CrossLink other)
		{
			if (ReferenceEquals(null, other)) return false;
			if (ReferenceEquals(this, other)) return true;
			return other.Mass.Equals(this.Mass) && Equals(other.ModType, this.ModType);
		}

		public override int GetHashCode()
		{
			unchecked
			{
				return (this.Mass.GetHashCode()*397) ^ this.ModType.GetHashCode();
			}
		}
	}
}
