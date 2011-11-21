using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CrossLinkingIMS.Control;

namespace CrossLinkingIMSConsole
{
	class ConsoleApplication
	{
		static void Main(string[] args)
		{
			CrossLinkingIMSController controller = new CrossLinkingIMSController();
			controller.Execute();
		}
	}
}
