#include "DivRateSlider.hpp"

DivRateSlider::DivRateSlider() :
		AbstractProposal(Category::DIVERGENCE_TIMES, "divRateSlider", 10, 1e-5,
				1e2, false), _savedContent
		{ }
{
}

DivRateSlider::~DivRateSlider()
{
}

double DivRateSlider::getNewProposal(double oldRate, Randomness &rand)
{
	/* the absolute value here ensures a positive value */
	return std::abs(rand.drawFromSlidingWindow(oldRate, 0.1));
}

void DivRateSlider::applyToState(TreeAln &traln, PriorBelief &prior,
		log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
	auto primParam = _allParams->at(_primParamIds[0]);
	_savedContent = primParam->extractParameter(traln);
	auto content = _savedContent;

	/* Select node at random */
	nat nodeToChange = rand.drawIntegerClosed(
			_savedContent.values.size() / 2 - 1);

	/* TODO-divtimes IGNORE ROOT NODES SO FAR! */
	if (traln.isRootChild(nodeToChange + 1))
		return;

	/* Propose a new value */
	double newRate = getNewProposal(_savedContent.values[nodeToChange], rand);

	/* Update branch length */
	double branchLength =
			content.branchLengths[nodeToChange].getLength().getValue();
	branchLength = exp(
			log(branchLength) * newRate / _savedContent.values[nodeToChange]);

	/* TODO-divtimes: discard if branch length is out of bounds! */
	if (branchLength
			&& (branchLength < BoundsChecker::zMin
					|| branchLength > BoundsChecker::zMax))
	{
		return;
	}

	content.values[nodeToChange] = newRate;
	content.branchLengths[nodeToChange].setLength(
			InternalBranchLength(branchLength));

#ifdef DIVTIME_DEBUG
	std::cout << "DIV-RATES-DBG: Propose new rate for " << nodeToChange << " --> " << _savedContent.values[nodeToChange] << " to " << newRate << std::endl;
	std::cout << "DIV-RATES-DBG: - - Branch changes from " << content.branchLengths[nodeToChange].getLength().getValue() << " to " << branchLength << std::endl;
#endif

	primParam->applyParameter(traln, content);
}

void DivRateSlider::evaluateProposal(LikelihoodEvaluator &evaluator,
		TreeAln &traln, const BranchPlain &branchSuggestion)
{
	auto prts = _allParams->at(_primParamIds[0])->getPartitions();
	evaluator.evaluatePartitionsWithRoot(traln, branchSuggestion, prts, true);
}

void DivRateSlider::resetState(TreeAln &traln)
{
#ifdef DIVTIME_DEBUG
	std::cout << "DIV-RATES-DBG: Reset rate proposal " << std::endl;
#endif
	_allParams->at(_primParamIds[0])->applyParameter(traln, _savedContent);
}

void DivRateSlider::autotune()
{
}

AbstractProposal* DivRateSlider::clone() const
{
	return new DivRateSlider(*this);
}

BranchPlain DivRateSlider::determinePrimeBranch(const TreeAln &traln,
		Randomness& rand) const
{
	return BranchPlain();
}

std::vector<nat> DivRateSlider::getInvalidatedNodes(const TreeAln &traln) const
{
	return
	{};
}

std::pair<BranchPlain, BranchPlain> DivRateSlider::prepareForSetExecution(
		TreeAln &traln, Randomness &rand)
{
	return std::make_pair(BranchPlain(0, 0), BranchPlain(0, 0));
}

void DivRateSlider::writeToCheckpointCore(std::ostream &out) const
{
}

void DivRateSlider::readFromCheckpointCore(std::istream &in)
{
}

