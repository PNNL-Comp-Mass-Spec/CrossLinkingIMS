using System.Collections.Generic;
using CrossLinkingIMS.Util;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS.Data
{
	// TODO: Explain the different ModTypes
	public enum ModType { None, Zero, One, Two, ZeroOne, ZeroTwo };

	/// <summary>
	/// Stores information about a single cross-link for either 1 or 2 peptides.
	/// </summary>
	public class CrossLink
	{
		public clsInSilicoDigest.PeptideInfoClass PeptideOne { get; private set; }
		public clsInSilicoDigest.PeptideInfoClass PeptideTwo { get; private set; }
		public double Mass { get; private set; }
		public ModType ModType { get; private set; }
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

		public override bool Equals(object obj)
		{
			if (ReferenceEquals(null, obj)) return false;
			if (ReferenceEquals(this, obj)) return true;
			if (obj.GetType() != typeof (CrossLink)) return false;
			return Equals((CrossLink) obj);
		}

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
