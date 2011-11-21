using System.Collections.Generic;
using System.Linq;

namespace CrossLinkingIMS.Data
{
	/// <summary>
	/// Class used for storing results for searching for the mass shifts of a CrossLink.
	/// </summary>
	public class MassShiftResult
	{
		/// <summary>
		/// This list stores a list of KeyValuePairs that contain the mass shift value searched for and if the mass value was found.
		/// </summary>
		public List<KeyValuePair<double, bool>> KvpList { get; set; }

		/// <summary>
		/// Default construtor.
		/// </summary>
		public MassShiftResult()
		{
			this.KvpList = new List<KeyValuePair<double, bool>>();
		}

		/// <summary>
		/// Calculates how many mass-shift values were found in the data for this specific CrossLinkResult object.
		/// </summary>
		/// <returns>The number of mass-shift values found in the data.</returns>
		public int CalculateMassShiftsFound()
		{
			return KvpList.Count(entry => entry.Value);
		}

		/// <summary>
		/// Determines whether the specified Object is equal to the current Object.
		/// The objects are equal if every entry in KvpList match.
		/// </summary>
		/// <param name="other">The specified MassShiftResult object.</param>
		/// <returns>true if the objects are equal, false otherwise</returns>
		public bool Equals(MassShiftResult other)
		{
			if (ReferenceEquals(null, other)) return false;
			if (ReferenceEquals(this, other)) return true;

			foreach (var keyValuePair in KvpList)
			{
				bool otherValue = other.KvpList.Where(x => x.Key == keyValuePair.Key).First().Value;
				if(otherValue != keyValuePair.Value)
				{
					return false;
				}
			}

			return true;
		}

		/// <summary>
		/// Determines whether the specified Object is equal to the current Object.
		/// The objects are equal if every entry in KvpList match.
		/// </summary>
		/// <param name="obj">The specified Object.</param>
		/// <returns>true if the objects are equal, false otherwise</returns>
		public override bool Equals(object obj)
		{
			if (ReferenceEquals(null, obj)) return false;
			if (ReferenceEquals(this, obj)) return true;
			if (obj.GetType() != typeof (MassShiftResult)) return false;
			return Equals((MassShiftResult) obj);
		}

		public override int GetHashCode()
		{
			unchecked
			{
				int hash = 17;

				foreach (var kvp in KvpList)
				{
					hash *= 23 + kvp.Key.GetHashCode() + kvp.Value.GetHashCode();
				}

				return hash;
			}
		}
	}
}
