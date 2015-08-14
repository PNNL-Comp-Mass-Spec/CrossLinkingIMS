using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using CrossLinkingIMS.Constants;
using CrossLinkingIMS.Data;
using ProteinDigestionSimulator;

namespace CrossLinkingIMS.Util
{
    /// <summary>
    /// Utility class containing function related to CrossLink.
    /// </summary>
    public class CrossLinkUtil
    {
        /// <summary>
        /// Fixed (static) mass difference
        /// </summary>
        public static double StaticDeltaMass { get; set; }

        /// <summary>
        /// When true, add a mass difference based on the number of carbon atoms 
        /// </summary>
        public static bool UseC13 
        {
            get;
            set;
        }

        /// <summary>
        /// When true, add a mass difference based on the number of nitrogen atoms 
        /// </summary>
        public static bool UseN15 
        {
            get;
            set;
        }

        private static readonly Dictionary<char, AminoAcid> m_aminoAcidDictionary;
        private static readonly char[] m_crossLinkCharacters = { 'K', 'S', 'T', 'Y' };

        /// <summary>
        /// Static constructor that will initialize the Amino Acid Dictionary
        /// </summary>
        static CrossLinkUtil()
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
        /// Generates a list of theoretical crosslinks based on the given protein sequence and list of peptides.
        /// </summary>
        /// <param name="peptideEnumerable">List of peptides to generate cross-links from.</param>
        /// <param name="proteinSequence">Protein sequence used.</param>
        /// <param name="proteinId">The identifier of the protein sequence used.</param>
        /// <returns>An IEnumerable of CrossLink objects.</returns>
        public static IEnumerable<CrossLink> GenerateTheoreticalCrossLinks(
          IEnumerable<clsInSilicoDigest.PeptideInfoClass> peptideEnumerable, 
          string proteinSequence, 
          string proteinId)
        {
            var crossLinkSet = new HashSet<CrossLink>();

            foreach (var firstPeptide in peptideEnumerable)
            {
                // First find all cross links using the first peptide + a null peptide
                foreach (var crossLink in FindCrossLinks(proteinId, firstPeptide, null, proteinSequence))
                {
                    crossLinkSet.Add(crossLink);
                }

                // Then find all cross links using the first peptide and every other peptide
                foreach (var secondPeptide in peptideEnumerable)
                {
                    foreach (var crossLink in FindCrossLinks(proteinId, firstPeptide, secondPeptide, proteinSequence))
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
        public static double CalculateMassShift(string sequenceString)
        {
            double numCarbon = 0;
            double numNitrogen = 0;

            foreach (var c in sequenceString)
            {
                var aminoAcid = m_aminoAcidDictionary[c];
                numCarbon += aminoAcid.NumCarbon;
                numNitrogen += aminoAcid.NumNitrogen;
            }

            double delMass = 0;

            if (UseC13)
                delMass += (numCarbon * (IsotopeConstants.CARBON_13_ISOTOPE - IsotopeConstants.CARBON_12_ISOTOPE));

            if (UseN15)
                delMass += (numNitrogen * (IsotopeConstants.NITROGEN_15_ISOTOPE - IsotopeConstants.NITROGEN_14_ISOTOPE));

            if (Math.Abs(StaticDeltaMass) > float.Epsilon)
                delMass += StaticDeltaMass;

            return delMass;

        }

        /// <summary>
        /// Finds all theoretical cross links given 2 peptides and a protein sequence.
        /// </summary>
        /// <param name="proteinId">The identifier of the protein used for cross linking.</param>
        /// <param name="firstPeptide">The first peptide used for cross linking.</param>
        /// <param name="secondPeptide">The second peptide used for cross linking. null if linking the first peptide to itself.</param>
        /// <param name="proteinSequence">Protein sequence used.</param>
        /// <returns>An IEnumerable of CrossLink objects.</returns>
        private static IEnumerable<CrossLink> FindCrossLinks(string proteinId, clsInSilicoDigest.PeptideInfoClass firstPeptide, clsInSilicoDigest.PeptideInfoClass secondPeptide, string proteinSequence)
        {
            var crossLinkList = new List<CrossLink>();

            // If 1 peptide
            if (secondPeptide == null)
            {
                // Create cross-link for unmodified peptide
                crossLinkList.Add(new CrossLink(proteinId, firstPeptide, null, firstPeptide.Mass, ModType.None));

                // Check for inter-linked peptides (will not always find)
                var peptideString = firstPeptide.SequenceOneLetter;
                if (peptideString.Last() == 'K')
                {
                    // Remove last character if it is a K
                    peptideString = peptideString.Substring(0, peptideString.Length - 1);
                }

                // Count the number of cross-link characters in the sequence
                var numCrossLinkCharacters = peptideString.Count(m_crossLinkCharacters.Contains);

                // If we are dealing with the peptide located at the very beginning of the protein sequence, then pretend we have an extra Lysine since we can cross-link with the first amino acid
                if (proteinSequence.StartsWith(peptideString))
                {
                    numCrossLinkCharacters++;
                }

                // If 0 Lysines are found, then we are done
                if (numCrossLinkCharacters == 0) return crossLinkList;

                // Type 0
                for (var i = 1; i <= numCrossLinkCharacters; i++)
                {
                    var modifiedMass = firstPeptide.Mass + (i * CrossLinkConstants.DEAD_END_MASS);
                    crossLinkList.Add(new CrossLink(proteinId, firstPeptide, null, modifiedMass, ModType.Zero));
                }

                // Type 1
                if (numCrossLinkCharacters >= 2)
                {
                    for (var i = 1; i <= numCrossLinkCharacters - 1; i++)
                    {
                        var modifiedMass = firstPeptide.Mass + (i * CrossLinkConstants.LINKER_MASS);
                        crossLinkList.Add(new CrossLink(proteinId, firstPeptide, null, modifiedMass, ModType.One));
                    }
                }

                // Type 0 and Type 1 mix
                if (numCrossLinkCharacters >= 3)
                {
                    for (var i = 1; i <= numCrossLinkCharacters / 2; i++)
                    {
                        var numLysinesLeft = numCrossLinkCharacters - (i * 2);
                        for (var j = 1; j <= numLysinesLeft; j++)
                        {
                            var modifiedMass = firstPeptide.Mass + (i * CrossLinkConstants.LINKER_MASS) + (j * CrossLinkConstants.DEAD_END_MASS);
                            crossLinkList.Add(new CrossLink(proteinId, firstPeptide, null, modifiedMass, ModType.ZeroOne));
                        }
                    }
                }
            }
            // If 2 peptides
            else
            {
                // First strip the last character from both peptide sequences if it is K, otherwise, leave it alone
                var firstPeptideString = firstPeptide.SequenceOneLetter;
                var secondPeptideString = secondPeptide.SequenceOneLetter;
                if (firstPeptideString.Last() == 'K')
                {
                    firstPeptideString = firstPeptideString.Substring(0, firstPeptideString.Length - 1);
                }
                if (secondPeptideString.Last() == 'K')
                {
                    secondPeptideString = secondPeptideString.Substring(0, secondPeptideString.Length - 1);
                }

                // Then count the number of Cross-Link characters in each sequence
                var firstPeptideNumCrossLinkCharacters = firstPeptideString.Count(m_crossLinkCharacters.Contains);
                var secondPeptideNumCrossLinkCharacters = secondPeptideString.Count(m_crossLinkCharacters.Contains);

                // If we are dealing with the peptide located at the very beginning of the protein sequence, then pretend we have an extra Lysine since we can cross-link with the first amino acid
                if (proteinSequence.StartsWith(firstPeptideString))
                {
                    firstPeptideNumCrossLinkCharacters++;
                }
                if (proteinSequence.StartsWith(secondPeptideString))
                {
                    secondPeptideNumCrossLinkCharacters++;
                }

                // If either peptide does not have a Lysine, then no cross-link is possible; exit
                if (firstPeptideNumCrossLinkCharacters == 0 || secondPeptideNumCrossLinkCharacters == 0)
                {
                    return crossLinkList;
                }

                // Add up the number of Lysines
                var numLysines = firstPeptideNumCrossLinkCharacters + secondPeptideNumCrossLinkCharacters;

                for (var i = 1; i <= numLysines / 2; i++)
                {
                    var numLysinesLeft = numLysines - (i * 2);
                    for (var j = 0; j <= numLysinesLeft; j++)
                    {
                        var modifiedMass = firstPeptide.Mass + secondPeptide.Mass + (i * CrossLinkConstants.LINKER_MASS) + (j * CrossLinkConstants.DEAD_END_MASS);

                        // Type 2
                        if (j == 0)
                        {
                            crossLinkList.Add(new CrossLink(proteinId, firstPeptide, secondPeptide, modifiedMass, ModType.Two));
                        }
                        // Type 2 and Type 0 mix
                        else
                        {
                            crossLinkList.Add(new CrossLink(proteinId, firstPeptide, secondPeptide, modifiedMass, ModType.ZeroTwo));
                        }
                    }
                }
            }

            return crossLinkList;
        }

        public static void OutputCrossLinks(IEnumerable<CrossLink> crossLinkEnumerable, FileInfo outputFileInfo)
        {
            TextWriter crossLinkWriter = new StreamWriter(outputFileInfo.FullName);
            crossLinkWriter.WriteLine("Index,Protein,Pep1,Pep2,ModType,TheoreticalMass");
            var index = 0;

            foreach (var crossLink in crossLinkEnumerable)
            {
                var outputLine = new StringBuilder();
                outputLine.Append(index++ + ",");
                outputLine.Append(crossLink.ProteinId + ",");
                outputLine.Append(crossLink.PeptideOne.SequenceOneLetter + ",");
                outputLine.Append((crossLink.PeptideTwo != null ? crossLink.PeptideTwo.SequenceOneLetter : "null") + ",");
                outputLine.Append(crossLink.ModType + ",");
                outputLine.Append(crossLink.Mass + "\n");

                crossLinkWriter.Write(outputLine);
            }

            crossLinkWriter.Close();
        }

        /// <summary>
        /// Writes the results of cross-link searching to a csv file.
        /// </summary>
        /// <param name="crossLinkResultEnumerable">The List of CrossLinkResults objects.</param>
        /// <param name="outputFileInfo"></param>
        public static void OutputCrossLinkResults(IEnumerable<CrossLinkResult> crossLinkResultEnumerable, FileInfo outputFileInfo)
        {
            TextWriter crossLinkWriter = new StreamWriter(outputFileInfo.FullName);
            crossLinkWriter.WriteLine("Index,Protein,Pep1,Pep2,ModType,TheoreticalMass,FeatureMass,PPMError,FeatureMz,ShiftedMassPep1,ShiftedMzPep1,ShiftedMassPep2,ShiftedMzPep2,ShiftedMassBoth,ShiftedMzBoth,ChargeState,LCScans,IMSScan,DriftTime,Abundance,FeatureIndex");
            var index = 0;

            var groupByCrossLinkQuery = crossLinkResultEnumerable.GroupBy(crossLinkResult => new
            {
                crossLinkResult.CrossLink,
                crossLinkResult.LcImsMsFeature.FeatureId,
                crossLinkResult.MassShiftResults
            })
                                                                 .Select(group => group.OrderBy(crossLinkResult => crossLinkResult.ScanLc));

            foreach (var group in groupByCrossLinkQuery)
            {
                var crossLinkResult = group.First();
                var crossLink = crossLinkResult.CrossLink;
                var feature = crossLinkResult.LcImsMsFeature;

                // Calculate PPM Error of the unshifted cross-link
                var ppmError = Math.Abs(crossLink.Mass - feature.MassMonoisotopic) / (crossLink.Mass / GeneralConstants.PPM_DIVISOR);

                var outputLine = new StringBuilder();
                outputLine.Append(index++ + ",");
                outputLine.Append(crossLink.ProteinId + ",");
                outputLine.Append(crossLink.PeptideOne.SequenceOneLetter + ",");
                outputLine.Append((crossLink.PeptideTwo != null ? crossLink.PeptideTwo.SequenceOneLetter : "null") + ",");
                outputLine.Append(crossLink.ModType + ",");
                outputLine.Append(crossLink.Mass + ",");
                outputLine.Append(feature.MassMonoisotopic + ",");
                outputLine.Append(ppmError + ",");
                outputLine.Append(feature.MzMonoisotopic + ",");

                // Iterate over the results for each mass-shift
                foreach (var massShiftResult in crossLinkResult.MassShiftResults.KvpList)
                {
                    // If the mass-shift value was found
                    if (massShiftResult.Value)
                    {
                        var shiftMass = massShiftResult.Key;
                        var shiftMz = (shiftMass / feature.ChargeState) + GeneralConstants.MASS_OF_PROTON;

                        // Output information about the mass-shift
                        outputLine.Append(shiftMass + "," + shiftMz + ",");
                    }
                    else
                    {
                        outputLine.Append("0,0,");
                    }

                    // If only 1 mass-shift, then the other 2 are N/A
                    if (crossLinkResult.MassShiftResults.KvpList.Count == 1)
                    {
                        outputLine.Append("N/A,N/A,N/A,N/A,");
                    }
                }

                outputLine.Append(feature.ChargeState + ",");

                foreach (var individualCrossLinkResult in group)
                {
                    outputLine.Append(individualCrossLinkResult.ScanLc + ";");
                }
                outputLine.Append(",");

                outputLine.Append(feature.ScanImsRep + ",");
                outputLine.Append(feature.DriftTime + ",");
                outputLine.Append(feature.Abundance + ",");
                outputLine.Append(feature.FeatureId + "\n");

                crossLinkWriter.Write(outputLine);
            }

            crossLinkWriter.Close();
        }
    }
}
