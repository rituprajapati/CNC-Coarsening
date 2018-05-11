This Implementation contains the idea of making CnC pure dataflow model (like TTG).
It is computing subtract(project(funcA), project(funcB)). In order to make CnC pure dataflow model,
for each of the structs Project and BinaryOp, we add:
  1) input_terminals (of type std::vector<CnC::item_collection<std::pair<int, int>, Node> *>) and
  2) outputTerminals (of type std::vector<OutputTerminal<std::pair<int, int>, Node>>)

OutputTerminal<K, V> structrue contains the following data members:
  1) CnC::item_collection<K, V> *item_collection, which is the output item_collection to put the output item.
  2) std::vector<CnC::tag_collection<K> *> next_op_tags, which is the tag_collection of the next operators.

By adding these two data members to the structs Project and BinaryOp and changing the constructor of the CnCContext struct, we have almost identical code in the execute() method of structs Project and BinaryOp (in file main.cc, which has the CnC implementation) with op() method defined in structs Project and BinaryOp (in file main_TTG.cc, which has the TTG implementation).

In order to compile this program on Seawulf, following commands need to be executed:
1) source <CNC_DIRECTORY>/bin/cncvars.sh
2) make

And to execute, the following command can be used:
./main MAX_LEVEL THRESHOLD

for example:
./main 5 0.000001

