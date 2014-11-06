""" General methods of general purpose
"""

def factorial(n):
    """ Returns factorial of n
    """
    if n == 0 or n == 1:
        return n
    else:
        return n*factorial(n-1)