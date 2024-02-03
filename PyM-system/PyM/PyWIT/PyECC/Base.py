# -*- coding: utf-8 -*-
"""
Created on Wed Jan 20 16:15:24 2016

@author: admin_local
"""

import PyM.PyWIT.PyECC.Base_type as Base_type

import copy

class Base(object):
    """Represents an abstract domain. """

    _group = False
    _ring = False
    _field = 0
    _mul_group = False
    _rel_structure = None
    _characteristic = None
    _poly_structure = False
    _matrix_structure = False
    _power_structure = False
    _fraction_structure = False
    _name = None
    _subDomain = None
    _base_Domain = None
    _ext_Domain = False
    _multi_variable = 0
    
    def subdomain(self):
        """returns its subdomain"""
        return self._subDomain
        
    def __rrshift__(self, a):
        raise NotImplementedError
    
    def __init__(self):
        raise NotImplementedError
    
    def _add_(self, a, b):
        """operates two elements using this domain operation +"""
        raise NotImplementedError
    
    def _sub_(self, a, b):
        """operates two elements using this domain operation -"""
        raise NotImplementedError
    
    def _mul_(self, a, b):
        """operates two elements using this domain operation *"""
        raise NotImplementedError
        
    def _inv_(self, a):
        """returns the inverse of the element a in this space"""
        raise NotImplementedError
        
    def _div_(self, a, b):
        """operates two elements using this domain operation /"""
        raise NotImplementedError
        
    def characteristic(self):
        """Return the characteristic of this domain. """
        raise NotImplementedError('characteristic()')
        
    def cardinal(self):
        """Return the cardinal of this domain. """
        raise NotImplementedError('cardinal()')
        
    def __eq__(self, other):
        raise NotImplementedError
        
    def __neq__(self, other):
        raise NotImplementedError
    
    def is_sub_element(self, a):
        """return if a is an element define in this space on any subspace"""
        if self.is_element(a):
            return True
        if self._subDomain.is_sub_element(a):
            return True
        return False
    
    def base_domain(self):
        """returns the first space of the chain of domains"""
        if self._base_Domain == None:
            return self
        return self._base_Domain
    
    def precedent_domain(self):
        """returns the last element of the chain of domains"""
        return self.subdomain()
        
    def is_sub_domain(self,other):
        """ returns if other is an element of the chain of domains"""
        if isinstance(other, Base):
            if self == other:
                return True
            return self.subdomain().is_sub_domain(other)
        return False
    
    def can_be_compared(self,other):
        if self.is_sub_domain(other):
            return True
        if other.is_sub_domain(self):
            return True
        A = self.K_()
        A2 = A.K_() 
        while(A!= A2):
            A = A2
            A2 = A.K_()
        B = other.K_()
        B2 = B.K_()
        while(B!= B2):
            B = B2
            B2 = B.K_()
        if A.is_sub_domain(B):
            return True
        if B.is_sub_domain(A):
            return True
        else:
            return False
    
    def union_domain(self,other): #Per ara no tenim en compte Q
        if self.is_sub_domain(other):
            return [self.cp(),0,0,0]
        if other.is_sub_domain(self):
            return [other.cp(),0,0,0]
        if self.can_be_compared(other):
            A = self
            KA = self.K_()
            lA = []
            B = other
            KB = other.K_()
            lB = []
            while (A!=KA):
                if A._poly_structure:
                    lA.append([1,A])
                elif A._power_structure:
                    lA.append([2,A,A._symbol])
                elif A._matrix_structure:
                    lA.append([3,A])
                else:
                    raise NotImplementedError
                A = A.subdomain()
            lA.reverse()
            while (B!=KB):
                if B._poly_structure:
                    lB.append([1,B])
                elif B._power_structure:
                    lB.append([2,B,B._symbol])
                elif B._matrix_structure:
                    lB.append([3,B])
                else:
                    raise NotImplementedError
                B = B.subdomain()
            lB.reverse()
            lC = []
            iA = 0
            iB = 0
            while(iA<len(lA) or iB<len(lB)):
                if iA == len(lA):
                    for i in range(iB,len(lB)):
                        lC.append([1]+lB[i])
                    iB = len(lB)
                    break
                if iB == len(lB):
                    for i in range(iA,len(lA)):
                        lC.append([0]+lA[i])
                    iA = len(lA)
                    break
                if lA[iA][0] < lB[iB][0]:
                    lC.append([0]+lA[iA])
                    iA += 1
                elif lB[iB][0] < lA[iA][0]:
                    lC.append([1]+lB[iB])
                    iB += 1
                else:
                    lC.append([2]+lA[iA])
                    iA += 1
                    iB += 1
            if KA.is_sub_domain(KB):
                D = KA
            else:
                D = KB
            symbUsed = []
            indsymbUsed = []
            for i in range(len(lC)):
                if lC[i][1] == 1:
                    D = lC[i][2]._create_Poly_Structure_(D)
                if lC[i][1] == 2:
                    if symbUsed.count(lC[i][3])>0:
                        s = symbUsed.index(lC[i][3])
                        symb = lC[i][3] + str(indsymbUsed[s])
                        indsymbUsed[s] += 1
                    else:
                        symb = lC[i][3]
                        symbUsed.append(lC[i][3])
                        indsymbUsed.append(1)
                    D = lC[i][2]._create_Power_Series_Structure_(D,symb)
                if lC[i][1] == 3:
                    D = lC[i][2]._create_Matrix_(D)
            return [D,lC,KA,KB]
        raise TypeError("These domains cannot be compared")
    
    def is_element(self, a):
        """returns if a is un element of THIS space"""
        if isinstance(a, Base_type.Base_type):
            if a._structure == self:
                return True
        return False
    
    def is_invertible(self, a):
        """returns if a is invertible in this space"""
        raise NotImplementedError
    
    def change_characteristic(self,n):
        """changes the characteristic of this space to n and transforms all objects that define this space"""
        raise NotImplementedError
    
    def transform_to_ZZ(self):
        """changes the characteristic of this space to 0 and transforms all objects that define this space"""
        raise NotImplementedError
    
    def set_name(self,name):
        """changes the name of the element"""
        self._name = name
    
    def get_name(self):
        """returns the name of the element"""
        return self._name
    
    def K_(self):
        raise NotImplementedError
    
    def cp(self):
        raise NotImplementedError 
    
    #Ficar les que surgeixin.