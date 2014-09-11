#include "DivergenceRates.hpp"

#include "Category.hpp"

DivergenceRates::DivergenceRates(nat id, nat idOfMyKind,
		std::vector<nat> partitions, nat numberOfTaxa) :
		AbstractParameter(Category::DIVERGENCE_RATES, id, idOfMyKind,
				partitions, 0), _rateAssignments(2 * numberOfTaxa - 2, 0), _rates(
				2 * numberOfTaxa - 2, 1.), _directedBranches(2 * numberOfTaxa - 2)
{
	/*
	 * we initialize as many different rate categories as branches in the tree, all to 1.0
	 */
	for (int i = 0; i < _rateAssignments.size(); i++)
	{
		_rateAssignments[i] = i;
		_rates[i] = 1;
	}
}

void DivergenceRates::initializeParameter(TreeAln& traln,
		const ParameterContent &content)
{
	assert(content.branchLengths.size() == traln.getNumberOfNodes());
	for (int i=0; i<traln.getNumberOfNodes();i++)
	{
		_directedBranches[i] = content.branchLengths[i];
	}

}

AbstractParameter* DivergenceRates::clone() const
{
	DivergenceRates * ndr = new DivergenceRates(*this);
	ndr->setRates(_rates);
	ndr->_rateAssignments = _rateAssignments;
	ndr->_directedBranches = _directedBranches;
	return ndr;
}

void DivergenceRates::applyParameter(TreeAln& traln,
		const ParameterContent &content)
{
	verifyContent(traln, content);

	for (int i = 0; i < _rates.size(); i++)
	{
		if (content.values[i] != _rates[i])
		{
			_rates[i] = content.values[i];
			traln.setBranch(content.branchLengths[i], this);
		}
	}

}

ParameterContent DivergenceRates::extractParameter(const TreeAln &traln) const
{
	auto content = ParameterContent();
	content.values = _rates;

	for (int i=0; i<_rates.size(); i++)
	{
		content.values.push_back(1.0 * _rateAssignments[i]);
		content.branchLengths.push_back(
					traln.getBranch(_directedBranches[i], this));
	}
	return content;
}

void DivergenceRates::printSample(std::ostream& fileHandle,
		const TreeAln &traln) const
{
	bool isFirst = true;
	for (auto v : _rates)
	{
		fileHandle << (isFirst ? "" : "\t") << v;
		isFirst = false;
	}
}

void DivergenceRates::printAllComponentNames(std::ostream &fileHandle,
		const TreeAln &traln) const
{
	bool isFirstG = true;
	for (nat i = 0; i < _rates.size(); ++i)
	{
		fileHandle << (isFirstG ? "" : "\t") << "r{";
		isFirstG = false;

		bool isFirst = true;
		for (auto &p : _partitions)
		{
			fileHandle << (isFirst ? "" : ",") << p;
			isFirst = false;
		}
		fileHandle << "}(" << i << ")";
	}
}

void DivergenceRates::verifyContent(const TreeAln &traln,
		const ParameterContent &content) const
{
	/*
	 * content must have
     * 1) rateValues (double) for every node
	 * 2) rateAssignments (nat) for every node
	 */
	assert(content.values.size() == 2*_rates.size());
	assert(content.branchLengths.size() == _directedBranches.size());
	for (int i=0;i<_directedBranches.size(); i++)
	{
		assert(content.branchLengths[i] == _directedBranches[i]);
	}
	for (int i=_rates.size(); i<2*_rates.size(); i++)
	{
		double intpart;
		assert(modf(content.values[i], &intpart) == 0.0);
	}
}

log_double DivergenceRates::getPriorValue(const TreeAln& traln) const
{

	/* TODO-divtimes check this */
	auto result = log_double::fromAbs(1);

	// that only partly makes sense ...
	// result *= _prior->getLogProb( _rates);

	return result;
}

void DivergenceRates::setPrior(const std::unique_ptr<AbstractPrior> &prior)
{
	// meh meh meh

	AbstractParameter::setPrior(prior);
	//_rates.resize(prior->getInitialValue().values.size());
}
