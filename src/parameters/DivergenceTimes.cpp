#include "DivergenceTimes.hpp"
#include "Category.hpp"

DivergenceTimes::DivergenceTimes(nat id, nat idOfMyKind, std::vector<nat> partitions, NodeAge age)
  : AbstractParameter(Category::DIVERGENCE_TIMES, id, idOfMyKind, partitions, 0 )
  , _nodeAge{age}
{
}

void DivergenceTimes::initializeParameter(TreeAln& traln,  const ParameterContent &content, bool root) {

	_nodeAge.setHeight(content.nodeAges[0].getHeight());
	_nodeAge.setPrimNode(content.nodeAges[0].getPrimNode());
	_nodeAge.setSecNode(content.nodeAges[0].getSecNode());

	if (root)
	{
		/* root case */
		if (traln.isTipBranch(_nodeAge)) {
			/* update the root branch */

		}
	}
	else
	{
		/* initialize parental branch length */
		if (!traln.isRootBranch(_nodeAge))
		{
			auto bl = traln.getBranch(_nodeAge, this);
			auto length = bl.getLength();

			double time = content.nodeAges[1].getHeight()
					- _nodeAge.getHeight();

			/* we assume that the current rate at initialization is 1.0 */
			double newLength = exp(-time / traln.getTrHandle().fracchange);

			auto newBranch = BranchLength(_nodeAge, newLength);
			traln.setBranch(newBranch, this);
		}

		/* initialize descendant branch length */
		auto descendants = traln.getDescendents(_nodeAge);
//		for (auto b : descendants)
//		{
//			if (traln.isTipBranch(b))
//			{
//				double time = _nodeAge.getHeight();
//				double newLength = exp(-time / traln.getTrHandle().fracchange);
//				auto newBranch = BranchLength(b, this);
//				newBranch.setLength(newLength);
//				traln.setBranch(newBranch, this);
//			}
//		}
		if (traln.isTipBranch(descendants.first))
		{
			double time = _nodeAge.getHeight();
			double newLength = exp(-time / traln.getTrHandle().fracchange);
			BranchLength newBranch = BranchLength(descendants.first, newLength);
			traln.setBranch(newBranch, this);
		}
		if (traln.isTipBranch(descendants.second))
		{
			double time = _nodeAge.getHeight();
			double newLength = exp(-time / traln.getTrHandle().fracchange);
			BranchLength newBranch = BranchLength(descendants.second,
					newLength);
			traln.setBranch(newBranch, this);
		}
	}
}

void DivergenceTimes::applyParameter(TreeAln& traln,  const ParameterContent &content)
{
}

 
ParameterContent DivergenceTimes::extractParameter(const TreeAln &traln)  const  
{
  auto content = ParameterContent();
  content.nodeAges.push_back(_nodeAge);
  return content;
}


void DivergenceTimes::printSample(std::ostream& fileHandle, const TreeAln &traln ) const 
{
		fileHandle << _nodeAge.getHeight();
}

 
void DivergenceTimes::printAllComponentNames(std::ostream &fileHandle, const TreeAln &traln) const  
{
		fileHandle << "t{";

		bool isFirst = true;
		for (auto &p : _partitions)
		{
			fileHandle << (isFirst ? "" : ",") << p;
			isFirst = false;
		}
		fileHandle << "}(" << _nodeAge << ")";
}

 
void DivergenceTimes::verifyContent(const TreeAln &traln, const ParameterContent &content) const
{
} 





log_double DivergenceTimes::getPriorValue(const TreeAln& traln) const
{
  // assert(0); return log_double::fromAbs(1); 
  // TODO should extract all branches and evaluate the prior ... not doing this here. .. 
  return log_double::fromAbs(1.);
} 
