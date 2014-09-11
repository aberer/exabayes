#include "DivergenceRates.hpp"

#include "Category.hpp"

DivergenceRates::DivergenceRates(nat id, nat idOfMyKind,
		std::vector<nat> partitions, nat numberOfTaxa) :
		AbstractParameter(Category::DIVERGENCE_RATES, id, idOfMyKind,
				partitions, 0), _rateAssignments(2 * numberOfTaxa - 2, 0), _rates(
				2 * numberOfTaxa - 2, 1.)
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

AbstractParameter* DivergenceRates::clone() const
{
	DivergenceRates * ndr = new DivergenceRates(*this);
	ndr->setRates(_rates);
	return ndr;
}

void DivergenceRates::applyParameter(TreeAln& traln,
		const ParameterContent &content)
{
	assert(content.values.size() == _rates.size());
	for (int i = 0; i < content.values.size(); i++)
	{
		if (content.nodeAges.size() == _rates.size())
		{
			std::cout << " RATEDBG " << content.values[i] << " " << content.nodeAges[i] << std::endl;
		}
		if (content.values[i] != _rates[i])
		{
			/* rate update */
			if (!traln.isRootBranch(i))
			{
				/* root branch */
			}
			else
			{
				/* TODO-divtimes We need to know which branch is pointing to the root! */
				//traln.getBranchesFromNode(i)
				//_rates[i] = content.values[i];
			}
			_rates[i] = content.values[i];
		}
	}
}

ParameterContent DivergenceRates::extractParameter(const TreeAln &traln) const
{
	auto content = ParameterContent();
	content.values = _rates;
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
//	assert(content.values.size() == traln.getNumberOfNodes());
	assert(content.values.size() == _rates.size());
}

log_double DivergenceRates::getPriorValue(const TreeAln& traln) const
{
	auto result = log_double::fromAbs(1);

	// that only partly makes sense ...
	// result *= _prior->getLogProb( _rates);

	return result;
}

void DivergenceRates::setPrior(const std::unique_ptr<AbstractPrior> &prior)
{
	// meh meh meh

	AbstractParameter::setPrior(prior);
	_rates.resize(prior->getInitialValue().values.size());
}
