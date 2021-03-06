#ifndef SRC_INPUT_STREAMINPUT_H
#define SRC_INPUT_STREAMINPUT_H

#include "ClusterEditingInput.h"
#include "ClusterEditingInstance.h"


namespace yskInput{

	/**
	 * Virtual class that serves as a template for any input format that is provided as an actual input stream (usually by means of an input file)
	 */
	class StreamInput : public ClusterEditingInput{
		public:

			StreamInput(ysk::ClusterEditingInstance* inst):ClusterEditingInput(inst){};

			/**
			 * Parses the input stream treating it as the input format represented by the respective class implementing this method
			 * @param is The input stream which is to be parsed
			 * @return True if the parsing was successful / False otherwise
			 */
			virtual bool parseInput(std::istream &is) = 0;

			enum Format{
				JENA = 0,
				SIF = 1,
				SIMROW = 2
			};


	};

}

#endif /* SRC_INPUT_STREAMINPUT_H */
