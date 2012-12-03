using System;
using System.Collections.Generic;
using System.Linq;
using CrossLinkingIMS.Constants;
using CrossLinkingIMS.Data;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS.Util
{
	/// <summary>
	/// Utility class for processing sequences.
	/// </summary>
	public class SequenceUtil
	{
		/// <summary>
		/// Given a protein sequence, a collection of peptides will be returned. 
		/// Assumes conventional trypsin digest and max 1 missed cleavage.
		/// </summary>
		/// <param name="proteinSequence">The protein sequence to digest.</param>
		/// <param name="digestionRule">Determines which digestion rule to use (fully tryptic, partially tryptic, no rules).</param>
		/// <param name="maxMissedCleavages">The maximum number of missed cleavages to consider.</param>
		/// <returns>An IEnumerable of Peptide objects.</returns>
		public static IEnumerable<clsInSilicoDigest.PeptideInfoClass> DigestProtein(string proteinSequence, clsInSilicoDigest.CleavageRuleConstants digestionRule, int maxMissedCleavages)
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
				CleavageRuleID = digestionRule,
				MaxMissedCleavages = maxMissedCleavages
			};

			parseProteinFile.DigestProteinSequence(proteinSequence, ref peptideArray, digestionOptions, "xlinkProt");
			
			return peptideArray;
		}
	}
}
