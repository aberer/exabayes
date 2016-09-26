#ifndef _PARAMETER_SET_HPP
#define _PARAMETER_SET_HPP

#include "AbstractParameter.hpp"
#include "BranchLengthsParameter.hpp"

/**
 *  @brief represents a set of parameters, that we integrate over
 *
 *  this object owns the respective parameters
 */
///////////////////////////////////////////////////////////////////////////////
//                               PARAMETER LIST                              //
///////////////////////////////////////////////////////////////////////////////
class ParameterList
{
    ///////////////////////////////////////////////////////////////////////////
    //                            PUBLIC INTERFACE                           //
    ///////////////////////////////////////////////////////////////////////////
public:
    // ________________________________________________________________________
    ParameterList(
        std::vector<std::unique_ptr<AbstractParameter> >params  =
            std::vector<std::unique_ptr<AbstractParameter> >{});
    // ________________________________________________________________________
    ~ParameterList(){}
    // ________________________________________________________________________
    ParameterList(
        ParameterList&& rhs) = default;
    // ________________________________________________________________________
    ParameterList(
        const ParameterList& rhs);
    // ________________________________________________________________________
    ParameterList&
                        operator
    =(
        ParameterList rhs);
    // ________________________________________________________________________
    friend void
                        swap(
        ParameterList& lhs,
        ParameterList& rhs);
    // ________________________________________________________________________
    // boring stuff to make it behave like a container
    const AbstractParameter*
                        at(
        int num) const {return _paramView.at(num); }
    // ________________________________________________________________________
    AbstractParameter*
                        at(
        int num){return _paramView.at(num); }
    // ________________________________________________________________________
    std::vector<AbstractParameter*>::const_iterator
                        begin()
    const
    {return _paramView.cbegin(); }
    // ________________________________________________________________________
    std::vector<AbstractParameter*>::const_iterator
                        end()
    const
    {return _paramView.cend(); }
    // ________________________________________________________________________
    std::vector<AbstractParameter*>::iterator
                        begin()
    {return _paramView.begin(); }
    // ________________________________________________________________________
    std::vector<AbstractParameter*>::iterator
                        end()
    {return _paramView.end(); }
    // ________________________________________________________________________
    AbstractParameter*
                        operator
    [](
        int num){return _paramView.at(num); }
    const AbstractParameter*
                        operator
    [](
        int num) const {return _paramView.at(num); }
    // ________________________________________________________________________
    friend std::ostream&
                        operator
    <<(
        std::ostream&        s,
        const ParameterList& c);
    // ________________________________________________________________________
    std::vector<AbstractParameter*>
                        getViewByCategory(
        Category cat) const;

    ///////////////////////////////////////////////////////////////////////////
    //                              PRIVATE DATA                             //
    ///////////////////////////////////////////////////////////////////////////
private:
    std::vector<std::unique_ptr<AbstractParameter> > _params;
    std::vector<AbstractParameter*>                  _paramView;
};


#endif
