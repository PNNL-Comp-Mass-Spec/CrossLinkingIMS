using System;
using System.Collections.Generic;
using System.Linq;
using CrossLinkingIMS.Constants;
using CrossLinkingIMS.Data;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS.Util
{
	/// <summary>
	/// Utility class for processing sequences and searching for cross-links.
	/// </summary>
	public class SequenceUtil
	{
		private static Dictionary<char, AminoAcid> m_aminoAcidDictionary;

		/// <summary>
		/// Static constructor that will initialize the Amino Acid Dictionary
		/// </summary>
		static SequenceUtil()
		{
			m_aminoAcidDictionary = new Dictionary<char, AminoAcid>
			{
				{'A', new AminoAcid {NumCarbon = 3, NumHydrogen = 5, NumNitrogen = 1, NumOxygen = 1}},
				{'C', new AminoAcid {NumCarbon = 3, NumHydrogen = 5, NumNitrogen = 1, NumOxygen = 1, NumSulfur = 1}},
				{'D', new AminoAcid {NumCarbon = 4, NumHydrogen = 5, NumNitrogen = 1, NumOxygen = 3}},
				{'E', new AminoAcid {NumCarbon = 5, NumHydrogen = 7, NumNitrogen = 1, NumOxygen = 3}},
				{'F', new AminoAcid {NumCarbon = 9, NumHydrogen = 9, NumNitrogen = 1, NumOxygen = 1}},
				{'G', new AminoAcid {NumCarbon = 2, NumHydrogen = 3, NumNitrogen = 1, NumOxygen = 1}},
				{'H', new AminoAcid {NumCarbon = 6, NumHydrogen = 7, NumNitrogen = 3, NumOxygen = 1}},
				{'I', new AminoAcid {NumCarbon = 6, NumHydrogen = 11, NumNitrogen = 1, NumOxygen = 1}},
				{'J', new AminoAcid {NumCarbon = 37, NumHydrogen = 36, NumNitrogen = 5, NumOxygen = 5, NumSulfur = 1, NumIron = 1}},
				{'K', new AminoAcid {NumCarbon = 6, NumHydrogen = 12, NumNitrogen = 2, NumOxygen = 1}},
				{'L', new AminoAcid {NumCarbon = 6, NumHydrogen = 11, NumNitrogen = 1, NumOxygen = 1}},
				{'M', new AminoAcid {NumCarbon = 5, NumHydrogen = 9, NumNitrogen = 1, NumOxygen = 1, NumSulfur = 1}},
				{'N', new AminoAcid {NumCarbon = 4, NumHydrogen = 6, NumNitrogen = 2, NumOxygen = 2}},
				{'P', new AminoAcid {NumCarbon = 5, NumHydrogen = 7, NumNitrogen = 1, NumOxygen = 1}},
				{'Q', new AminoAcid {NumCarbon = 5, NumHydrogen = 8, NumNitrogen = 2, NumOxygen = 2}},
				{'R', new AminoAcid {NumCarbon = 6, NumHydrogen = 12, NumNitrogen = 4, NumOxygen = 1}},
				{'S', new AminoAcid {NumCarbon = 3, NumHydrogen = 5, NumNitrogen = 1, NumOxygen = 2}},
				{'T', new AminoAcid {NumCarbon = 4, NumHydrogen = 7, NumNitrogen = 1, NumOxygen = 2}},
				{'V', new AminoAcid {NumCarbon = 5, NumHydrogen = 9, NumNitrogen = 1, NumOxygen = 1}},
				{'W', new AminoAcid {NumCarbon = 11, NumHydrogen = 10, NumNitrogen = 2, NumOxygen = 1}},
				{'Y', new AminoAcid {NumCarbon = 9, NumHydrogen = 9, NumNitrogen = 1, NumOxygen = 2}}
			};
		}

		/// <summary>
		/// Given a protein sequence, a collection of peptides will be returned. 
		/// Assumes conventional trypsin digest and max 1 missed cleavage.
		/// </summary>
		/// <param name="proteinSequence">The protein sequence to digest.</param>
		/// <returns>An IEnumerable of Peptide objects.</returns>
		public static IEnumerable<clsInSilicoDigest.PeptideInfoClass> DigestProtein(string proteinSequence)
		{
			clsParseProteinFile parseProteinFile = new clsParseProteinFile
			{
				AssumeFastaFile = true,
				AssumeDelimitedFile = false,
				ComputeProteinMass = true,
				CreateProteinOutputFile = false,
				CreateDigestedProteinOutputFile = false,
				GenerateUniqueIDValuesForPeptides = true
			};

			clsInSilicoDigest.PeptideInfoClass[] peptideArray = new clsInSilicoDigest.PeptideInfoClass[1];

			clsInSilicoDigest.DigestionOptionsClass digestionOptions = new clsInSilicoDigest.DigestionOptionsClass
			{
				CleavageRuleID = clsInSilicoDigest.CleavageRuleConstants.ConventionalTrypsin,
				MaxMissedCleavages = 1
			};

			parseProteinFile.DigestProteinSequence(proteinSequence, ref peptideArray, digestionOptions, "xlinkProt");

			foreach (clsInSilicoDigest.PeptideInfoClass peptide in peptideArray)
			{
				Console.WriteLine(peptide.SequenceOneLetter);
			}

			return peptideArray;
		}

		/// <summary>
		/// Generates a list of theoretical crosslinks based on the given protein sequence and list of peptides.
		/// </summary>
		/// <param name="peptideEnumerable">List of peptides to generate cross-links from.</param>
		/// <param name="proteinSequence">Protein sequence used.</param>
		/// <returns>An IEnumerable of CrossLink objects.</returns>
		public static IEnumerable<CrossLink> GenerateTheoreticalCrossLinks(IEnumerable<clsInSilicoDigest.PeptideInfoClass> peptideEnumerable, String proteinSequence)
		{
			HashSet<CrossLink> crossLinkSet = new HashSet<CrossLink>();

			foreach (clsInSilicoDigest.PeptideInfoClass firstPeptide in peptideEnumerable)
			{
				// First find all cross links using the first peptide + a null peptide
				foreach (CrossLink crossLink in FindCrossLinks(firstPeptide, null, proteinSequence))
				{
					crossLinkSet.Add(crossLink);
				}

				// Then find all cross links using the first peptide and every other peptide
				foreach (clsInSilicoDigest.PeptideInfoClass secondPeptide in peptideEnumerable)
				{
					foreach (CrossLink crossLink in FindCrossLinks(firstPeptide, secondPeptide, proteinSequence))
					{
						crossLinkSet.Add(crossLink);
					}
				}
			}

			return crossLinkSet;
		}

		/// <summary>
		/// Given a sequence of amino acids, calculates the mass shift.
		/// </summary>
		/// <param name="sequenceString">The sequence used for calculating the mass shift.</param>
		/// <returns>The calculated mass shift, in daltons.</returns>
		public static double CalculateMassShift(String sequenceString)
		{
			double numCarbon = 0;
			double numNitrogen = 0;

			foreach (char c in sequenceString)
			{
				AminoAcid aminoAcid = m_aminoAcidDictionary[c];
				numCarbon += aminoAcid.NumCarbon;
				numNitrogen += aminoAcid.NumNitrogen;
			}

			return (numCarbon * (13.00335484 - 12.0)) + (numNitrogen * (15.00010897 - 14.00307401));
		}

		/// <summary>
		/// Finds all theoretical cross links given 2 peptides and a protein sequence.
		/// </summary>
		/// <param name="firstPeptide">The first peptide used for cross linking.</param>
		/// <param name="secondPeptide">The second peptide used for cross linking. null if linking the first peptide to itself.</param>
		/// <param name="proteinSequence">Protein sequence used.</param>
		/// <returns>An IEnumerable of CrossLink objects.</returns>
		private static IEnumerable<CrossLink> FindCrossLinks(clsInSilicoDigest.PeptideInfoClass firstPeptide, clsInSilicoDigest.PeptideInfoClass secondPeptide, String proteinSequence)
		{
			List<CrossLink> crossLinkList = new List<CrossLink>();

			// If 1 peptide
			if(secondPeptide == null)
			{
				// Create cross-link for unmodified peptide
				crossLinkList.Add(new CrossLink(firstPeptide, null, firstPeptide.Mass, ModType.None));

				// Check for inter-linked peptides (will not always find)
				string peptideString = firstPeptide.SequenceOneLetter.Substring(0, firstPeptide.SequenceOneLetter.Length - 1); // Remove last character
				int numLysines = peptideString.Count(c => c == 'K');

				// If we are dealing with the peptide located at the very beginning of the protein sequence, then pretend we have an extra Lysine since we can cross-link with the first amino acid
				if(proteinSequence.StartsWith(peptideString))
				{
					numLysines++;
				}

				// If 0 Lysines are found, then we are done
				if (numLysines == 0) return crossLinkList;

				// Type 0
				for(int i = 1; i <= numLysines; i++)
				{
					double modifiedMass = firstPeptide.Mass + (i * CrossLinkConstants.DEAD_END_MASS);
					crossLinkList.Add(new CrossLink(firstPeptide, null, modifiedMass, ModType.Zero));
				}

				// Type 1
				if(numLysines >= 2)
				{
					for(int i = 1; i <= numLysines - 1; i++)
					{
						double modifiedMass = firstPeptide.Mass + (i * CrossLinkConstants.LINKER_MASS);
						crossLinkList.Add(new CrossLink(firstPeptide, null, modifiedMass, ModType.One));
					}
				}

				// Type 0 and Type 1 mix
				if(numLysines >= 3)
				{
					for(int i = 1; i <= numLysines/2; i++)
					{
						int numLysinesLeft = numLysines - (i * 2);
						for(int j = 1; j <= numLysinesLeft; j++)
						{
							double modifiedMass = firstPeptide.Mass + (i * CrossLinkConstants.LINKER_MASS) + (j * CrossLinkConstants.DEAD_END_MASS);
							crossLinkList.Add(new CrossLink(firstPeptide, null, modifiedMass, ModType.ZeroOne));
						}
					}
				}
			}
			// If 2 peptides
			else
			{
				// First strip the last character from both peptide sequences
				string firstPeptideString = firstPeptide.SequenceOneLetter.Substring(0, firstPeptide.SequenceOneLetter.Length - 1);
				string secondPeptideString = secondPeptide.SequenceOneLetter.Substring(0, secondPeptide.SequenceOneLetter.Length - 1);

				// Then count the number of Lysines
				int numLysines = 0;
				numLysines += firstPeptideString.Count(c => c == 'K');
				numLysines += secondPeptideString.Count(c => c == 'K');

				// If we are dealing with the peptide located at the very beginning of the protein sequence, then pretend we have an extra Lysine since we can cross-link with the first amino acid
				if (proteinSequence.StartsWith(firstPeptideString))
				{
					numLysines++;
				}
				if (proteinSequence.StartsWith(secondPeptideString))
				{
					numLysines++;
				}

				for (int i = 1; i <= numLysines / 2; i++)
				{
					int numLysinesLeft = numLysines - (i * 2);
					for (int j = 0; j <= numLysinesLeft; j++)
					{
						double modifiedMass = firstPeptide.Mass + secondPeptide.Mass + (i * CrossLinkConstants.LINKER_MASS) + (j * CrossLinkConstants.DEAD_END_MASS);

						if(j == 0)
						{
							crossLinkList.Add(new CrossLink(firstPeptide, secondPeptide, modifiedMass, ModType.Two));
						}
						else
						{
							crossLinkList.Add(new CrossLink(firstPeptide, secondPeptide, modifiedMass, ModType.ZeroTwo));
						}
					}
				}
			}

			return crossLinkList;
		}

		
	}
}
