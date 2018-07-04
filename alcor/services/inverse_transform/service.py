import operator
from typing import Iterator

import numpy as np
import sympy as sym
from scipy.special import lambertw
from sympy.functions.elementary.piecewise import ExprCondPair


def sample(pdf: sym.Function,
           *,
           size: int) -> np.array:
    """
    Generates random values following the given distribution
    :param pdf: input Probability Density Function (PDF)
    :param size: number of generated values
    """
    if not isinstance(pdf, sym.Piecewise):
        raise ValueError("PDF must be constructed by sympy.Piecewise")

    pdf_functions = map(operator.attrgetter('func'),
                        pdf.atoms(sym.Function))
    if sym.re in pdf_functions:
        error_message = ("Using sympy.Abs or sympy.re is not supported "
                         "due to not implemented computing of their integrals "
                         "within SymPy. Split the relevant expression.")
        raise NotImplementedError(error_message)

    # The following is used in order to prevent an error
    # when using PDF in a form of, for example, x**-2.5.
    # For more details see:
    # https://stackoverflow.com/questions/50543587/integrating-piecewise-with-irrational-exponent-gives-error
    pdf = sym.nsimplify(pdf)

    x = pdf.free_symbols.pop()
    y = sym.Dummy('y')

    cdf = sym.integrate(pdf, (x, -sym.oo, y))
    # The following is used in order to prevent
    # long erroneous polynomials
    # when calculating PDF in a form of, for example,  x**-2.5
    # Beware that this will add too much precision. Bug.
    # Issue submitted: https://github.com/sympy/sympy/issues/14787
    cdf = cdf.evalf()

    eq = sym.Eq(x, cdf)

    # TODO: Use solveset when it will be able to deal with LambertW
    # With default rational == True, there will be an error
    # as 'solve' doesn't play along with Piecewise.
    # Related issue: https://github.com/sympy/sympy/issues/12024
    inverse_solutions = sym.solve(eq, y, rational=False)
    # Sometimes, especially for exponents,
    # there are garbage solutions with imaginary parts:
    # https://github.com/sympy/sympy/issues/9973
    inverse_solutions = filter(is_real, inverse_solutions)

    # As, for some reason, 'solve' returns a list of Piecewise's,
    # it's necessary to collect them back together.
    # Related issue: https://github.com/sympy/sympy/issues/14733
    inverse_cdf = recreate_piecewise(inverse_solutions)
    # If inverse CDF will contain LambertW function,
    # we must change its branch. For more details, see:
    # https://stackoverflow.com/questions/49817984/sympy-solve-doesnt-give-one-of-the-solutions-with-lambertw
    functions = map(operator.attrgetter('func'),
                    inverse_cdf.atoms(sym.Function))
    if sym.LambertW in functions:
        inverse_cdf = replace_lambertw_branch(inverse_cdf)
        # This is to prevent LambertW giving ComplexWarning after lambdifying
        inverse_cdf = sym.re(inverse_cdf)

    max_value = cdf.args[-1][0]

    # Warnings can happen with exponents in PDF:
    # https://github.com/sympy/sympy/issues/14789
    lambda_function = sym.lambdify(args=x,
                                   expr=inverse_cdf,
                                   modules=[{'LambertW': lambertw}, 'numpy'])
    return lambda_function(np.random.uniform(high=max_value,
                                             size=size))


def is_real(expression: sym.Expr) -> bool:
    """Checks if expression doesn't contain imaginary part with sympy.I"""
    return sym.I not in expression.atoms()


def recreate_piecewise(functions: Iterator[ExprCondPair]) -> sym.Piecewise:
    """
    Collects Piecewise from list of unsorted Piecewise's,
    ignoring parts with NaNs.
    Solution for the issue: https://github.com/sympy/sympy/issues/14733
    See also question on SO:
    https://stackoverflow.com/questions/50428912/how-to-get-sorted-exprcondpairs-in-a-piecewise-function-that-was-obtained-from
    """
    def remove_nans(expression_condition: ExprCondPair) -> ExprCondPair:
        return expression_condition.args[0]

    def right_hand_number(solution: ExprCondPair) -> sym.S:
        return solution[1].args[1]

    solutions = sorted(map(remove_nans, functions),
                       key=right_hand_number)
    return sym.Piecewise(*solutions)


def to_lower_lambertw_branch(*args) -> sym.Function:
    """
    Wraps the first argument from a given list of arguments
    as a lower branch of LambertW function.
    :return: lower LambertW branch
    """
    return sym.LambertW(args[0], -1)


def replace_lambertw_branch(expression: sym.Expr) -> sym.Expr:
    """
    Replaces upper branch of LambertW function with the lower one.
    For details of the bug see:
    https://stackoverflow.com/questions/49817984/sympy-solve-doesnt-give-one-of-the-solutions-with-lambertw
    Solution is based on the 2nd example from:
    http://docs.sympy.org/latest/modules/core.html?highlight=replace#sympy.core.basic.Basic.replace
    :return: expression with replaced LambertW branch by a lower one
    """
    return expression.replace(sym.LambertW,
                              to_lower_lambertw_branch)
