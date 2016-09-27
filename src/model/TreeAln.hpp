/**
 * @brief represents a tree and the associated alignment
 *
 */

#ifndef _TREEALN_H
#define _TREEALN_H

#include "extensions.hpp"
#include "ArrayReservoir.hpp"
#include "GlobalVariables.hpp"

#include "ProtModel.hpp"
#include "BranchPlain.hpp"
#include "FlagType.hpp"

#include "BranchLengthResource.hpp"
#include "Partition.hpp"

#include <vector>
#include <iostream>
#include <array>
#include <memory>

class BranchLength;
class BranchLengths;
class AbstractPrior;
class Randomness;
class TreePrinter;
class AbstractParameter;
class ByteFile;


///////////////////////////////////////////////////////////////////////////////
//                                  TREE ALN                                 //
///////////////////////////////////////////////////////////////////////////////
class TreeAln
{
    friend class BranchLengthResource;

    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    TreeAln(
        size_t numTax,
        bool   isSaveMemorySEV);
    // ________________________________________________________________________
    TreeAln(
        const TreeAln& rhs);
    // ________________________________________________________________________
    TreeAln(
        TreeAln&& rhs) = default;
    // ________________________________________________________________________
    ~TreeAln(){}
    // ________________________________________________________________________
    TreeAln&                                          operator=(
        TreeAln rhs);
    // ________________________________________________________________________
    friend void                                       swap(
        TreeAln& lhs,
        TreeAln& rhs);
    // ________________________________________________________________________
    friend std::ostream&                              operator<<(
        std::ostream&   out,
        const TreeAln&  traln);
    // ________________________________________________________________________
    Partition&                                        getPartition(
        nat model);
    // ________________________________________________________________________
    const Partition&                                  getPartition(
        nat model)  const;
    // ________________________________________________________________________
    /**
     *  @brief get the internal raxml tree representation
     */
    pllInstance&                                      getTrHandle()
    {return _tr; }
    // ________________________________________________________________________
    const pllInstance&                                getTrHandle() const
    {return _tr; }
    // ________________________________________________________________________
    /**
     *  @brief frees all likelihood arrays
     */
    void                                              clearMemory(
        ArrayReservoir&arrayReservoir);
    // ________________________________________________________________________
    void                                              clearMemory();
    // ________________________________________________________________________
    /**
     *  @brief get the number of branches in the tree (not counting
     * per-partition branch lengths)
     */
    nat                                               getNumberOfBranches()
    const
    {return getNumberOfNodes() - 1; }
    // ________________________________________________________________________
    /**
     *  @brief get the number of partitions in the alignment
     */
    size_t                                            getNumberOfPartitions()
    const;
    // ________________________________________________________________________
    /**
     *  @brief get the number of taxa in the tree
     */
    nat                                               getNumberOfTaxa()
    const {return getTrHandle().mxtips; }
    // ________________________________________________________________________
    /**
     *  @brief get the number of nodes in the unrooted tree
     */

    nat                                               getNumberOfNodes()
    const
    {nat numTax = getNumberOfTaxa();  return 2 * numTax - 2;  }
    // excluding the virtual root
    // ________________________________________________________________________
    /**
     *  @brief get the substitution matrix for partition "model"
     */
    std::vector<double>                               getRevMat(
        nat  model,
        bool isRaw) const;
    // ________________________________________________________________________
    /**
     *  @brief gets the state frequencies for partition "model"
     */
    std::vector<double>                               getFrequencies(
        nat model) const;
    // ________________________________________________________________________
    /**
     * @brief gets the alpha parameter of the gamma distribution
     */
    double                                            getAlpha(
        nat model) const
    {
        auto&p = getPartition(model);
        return p.getAlpha();
    }
    // ________________________________________________________________________
    BranchPlain                                       getThirdBranch(
        const BranchPlain& oneBranch,
        const BranchPlain& otherBranch) const;
    // ________________________________________________________________________
    // not so happy with that...
    log_double                                        getLikelihood()
    const
    {
        return log_double::fromLog(getTrHandle().likelihood);
    }
    // ________________________________________________________________________
    void                                              setLikelihood(
        log_double val)
    {
        getTrHandle().likelihood = val.getRawLog();
    }
    // ________________________________________________________________________
    /**
     *  @brief indicates whether a nodepointer is a tip
     */
    bool                                              isTipNode(
        nodeptr p) const {return p->number <=  getTrHandle().mxtips; }
    // ________________________________________________________________________
    bool                                              isTipNode(
        nat aNode) const {return int(aNode) <=  getTrHandle().mxtips; }
    // ________________________________________________________________________
    /**
     *  @brief gets the branch from a node pointer (including branch length)
     *
     *  @notice if only one parameter is provided, the resulting branch
     *  will ONLY contain the branch length for this one parameter
     *  (bogus parameters for other branch lengths, if you have
     *  per-partition branch lengths). For setting this branch, the same
     *  parameter must be used (and only this one).
     *
     *  If you call this function with an existing Branch b, this
     *  function is usefull to get the actual branch length.
     */
    BranchLength                                      getBranch(
        const BranchPlain&      branch,
        const AbstractParameter*param) const;
    // ________________________________________________________________________
    BranchLengths                                     getBranch(
        const BranchPlain&                    branch,
        const std::vector<AbstractParameter*>&params) const;
    // ________________________________________________________________________
    BranchLength                                      getBranch(
        const nodeptr&          branch,
        const AbstractParameter*param) const;
    // ________________________________________________________________________
    BranchLengths                                     getBranch(
        const nodeptr&                        branch,
        const std::vector<AbstractParameter*>&params) const;
    // ________________________________________________________________________
    // improved interface for topological moves
    void                                              insertSubtree(
        const BranchPlain&                    subtree,
        const BranchPlain&                    insertBranch,
        const BranchLengths&                  branchToCreate,
        const std::vector<AbstractParameter*>&params);
    // ________________________________________________________________________
    auto                                              pruneSubtree(
        const BranchPlain&                    subtree,
        const BranchPlain&                    prunedBranch,
        const std::vector<AbstractParameter*>&params)
        ->std::tuple<BranchLengths, BranchPlain>;
    // ________________________________________________________________________
    /**
     *  @brief gets a nodepointer with specified id
     */
    nodeptr                                           getNode(
        nat elem) const;
    // ________________________________________________________________________
    std::vector<BranchPlain>                          getBranchesFromNode(
        nat aNode) const;
    // ________________________________________________________________________
    /**
     *  @brief extract all branches from the tree (including branch lengths)
     */
    // template<typename RESULT>
    std::vector<BranchLength>                         extractBranches(
        const AbstractParameter* param) const;
    // ________________________________________________________________________
    std::vector<BranchLengths>                        extractBranches(
        const std::vector<AbstractParameter*>&param) const;
    // ________________________________________________________________________
    std::vector<BranchPlain>                          extractBranches() const;
    /**
     *  @brief gets the number of inner nodes in the tree
     */
    // ________________________________________________________________________
    nat                                               getNumberOfInnerNodes(
        bool rooted) const;
    // ________________________________________________________________________
    /**
     *  @brief gets the mean substitution rate overall specified partitions
     */
    // double getMeanSubstitutionRate(const std::vector<nat> &partitions) const
    // ;
    /**
     *  @brief unlinks a node
     */
    void                                              unlinkNode(
        nodeptr p);
    // ________________________________________________________________________
    /**
     *  @brief gets the three nodes adjacent to the given node
     */
    std::vector<nat>                                  getNeighborsOfNode(
        nat node) const;
    // ________________________________________________________________________
    /**
     *  @brief prunes the node from the tree
     */
    void                                              detachNode(
        nodeptr p);
    // ________________________________________________________________________
    partitionList&                                    getPartitionsHandle();
    // ________________________________________________________________________
    const partitionList&                              getPartitionsHandle()
    const;
    // ________________________________________________________________________
    void                                              setPartitions(
        const std::vector<Partition>&p,
        bool                         initial);
    // ________________________________________________________________________
    BranchPlain                                       getAnyBranch()
    const;
    // ________________________________________________________________________
    /**
     * @brief sets the frequencies. Format is important, frequencies must add
     * up to 1.0
     */
    void                                              setFrequencies(
        const std::vector<double>&values,
        nat                       model);
    // ________________________________________________________________________
    /**
     *  @brief sets the parameters. Format is important, last rate must be 1.0
     */
    void                                              setRevMat(
        const std::vector<double>&values,
        nat                       model,
        bool                      isRaw);
    // ________________________________________________________________________
    /**
     *  @brief sets the alpha for partition "model"
     */
    void                                              setAlpha(
        double alpha,
        nat    model);
    // ________________________________________________________________________
    /**
     *  @brief sets a branch. Topology is NOT modified!
     */
    void                                              setBranch(
        const BranchLengths&                 b,
        const std::vector<AbstractParameter*>params);
    // ________________________________________________________________________
    void                                              setBranch(
        const BranchLength&      branch,
        const AbstractParameter* param);
    // ________________________________________________________________________
    /**
     *  @brief hooks up two nodes. Branch lengths must be set
     *  separately.
     */
    void                                              clipNode(
        nodeptr p,
        nodeptr q,
        double* z = nullptr);
    // ________________________________________________________________________
    bool                                              exists(
        const BranchPlain&branch) const;
    // ________________________________________________________________________
    /**
     *  @brief hooks up two nodes with default branch length
     */
    void                                              clipNodeDefault(
        nodeptr p,
        nodeptr q);
    // ________________________________________________________________________
    /**
     * @brief resets/destroys the topology in the tree
     */
    void                                              unlinkTree();
    // ________________________________________________________________________
    /**
     *  @brief gets the maximum length of paths below this branch
     */
    nat                                               getDepth(
        const BranchPlain&b) const;
    // ________________________________________________________________________
    /**
     *  @brief gets the longest path
     */
    std::vector<nat>                                  getLongestPathBelowBranch(
        const BranchPlain&b) const;
    // ________________________________________________________________________
    bool                                              isTipBranch(
        const BranchPlain&branch) const;
    // ________________________________________________________________________
    nodeptr                                           findNodePtr(
        const BranchPlain&branch) const;
    // ________________________________________________________________________
    /**
     * @brief gets a node with given id that is not connected to the tree right
     * now
     */
    nodeptr                                           getUnhookedNode(
        int number);
    // ________________________________________________________________________
    /**
     *  @brief gets the branches that are below one branch
     *
     *  This function is handy for traversing the tree and relying less
     *  on raw pointers. For traversing, it is often necessary to invert
     *  the resulting branch.
     */
    std::pair<BranchPlain,
              BranchPlain>                            getDescendents(
        const BranchPlain&b) const;
    // ________________________________________________________________________
    std::vector<bool>                                 getExecModel() const;
    // ________________________________________________________________________
    std::vector<log_double>                           getPartitionLnls() const;
    // ________________________________________________________________________
    void                                              setPartitionLnls(
        const std::vector<log_double>partitionLnls);
    // ________________________________________________________________________
    void                                              setExecModel(
        const std::vector<bool>&modelInfo);
    // ________________________________________________________________________
    nat                                               getNumberOfAssignedSites(
        nat model) const;
    // ________________________________________________________________________
    const std::vector<std::string>&                   getTaxa() const
    {return _taxa; }
    // ________________________________________________________________________
    void                                              setTaxa(
        std::vector<std::string>taxa){_taxa = taxa; }
    // ________________________________________________________________________
    std::vector<BranchPlain>                          getBranchesByDistance(
        const BranchPlain& branch,
        nat                distance,
        bool               bothSides) const;
    // ________________________________________________________________________
    void                                              setProteinModel(
        int       part,
        ProtModel model);
    // ________________________________________________________________________
    ProtModel                                         getProteinModel(
        int part) const;
    // ________________________________________________________________________
    void                                              setBranchUnchecked(
        const BranchLength&bl);
    // ________________________________________________________________________
    void                                              copyTopologyAndBl(
        const TreeAln&rhs);
    // ________________________________________________________________________
    void                                              setBranchLengthResource(
        BranchLengthResource bls);
    // ________________________________________________________________________
    bool                                              isSaveMemorySEV() const
    {return _isSaveMemorySEV; }
    // ________________________________________________________________________
    void                                              createCaterpillar();
    // ________________________________________________________________________
    void                                              initRevMat(
        nat model);             // these functions are not needed any more:
                                // directly use the respective setter function
    // ________________________________________________________________________
    void                                              discretizeGamma(
        nat model);
    // ________________________________________________________________________
    void                                              initialize(
        size_t numTax);
    // ________________________________________________________________________
    nat                                               addNodeToPartialTree(
        nat id,
        nat curRoot,
        nat outerCtr);
    // ________________________________________________________________________
    bool                                              isRooted(
        void) const;
    // ________________________________________________________________________
    void                                              setRootBranch(
        const BranchPlain&rb);
    // ________________________________________________________________________
    BranchPlain                                       getRootBranch() const;
    // ________________________________________________________________________
    bool                                              isRootChild(
        const nat nodeId) const;
    // ________________________________________________________________________
    bool                                              isRootBranch(
        const BranchPlain&rb) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    // ATTRIBUTES
    std::vector<Partition>     _partitions;
    partitionList              _partitionListResource;
    std::vector<pInfo*>        _partitionPtrs;

    pllInstance                _tr;
    std::vector<std::string>   _taxa;
    BranchLengthResource       _bls;

    // tree resources
    std::vector<traversalInfo> _traversalInfo;
    std::vector<node>          _nodes;
    std::vector<nodeptr>       _nodeptrs;

    std::vector<int>           _ti;
    std::vector<boolean>       _execModel; // for the traveral descriptor
    std::vector<nat>           _parsimonyScore;

    bool                       _isSaveMemorySEV;
    // in case we are working with a rooted tree
    BranchPlain                _root;
};


#endif
