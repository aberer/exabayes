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

void DivRateSlider::applyToState(TreeAln &traln, PriorBelief &prior,
		log_double &hastings, Randomness &rand, LikelihoodEvaluator& eval)
{
	auto primParam = _allParams->at(_primParamIds[0]);
	_savedContent = primParam->extractParameter(traln);
	//auto content = ParameterContent();
	auto content = _savedContent;
	content.nodeAges.resize(_secParamIds.size());
	for (int id : _secParamIds)
	{
		auto secParam = _allParams->at(_secParamIds[id]);
		content.nodeAges[secParam->getIdOfMyKind()] =
				secParam->extractParameter(traln).nodeAges[0];
	}
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

