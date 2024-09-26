
def recognize_impurity(name): 
    "Recognize the impurity based on a string"
    
    # Initiate the impurity
    impurity = None
    
    # Recognize the impurity from a string 
    if name=="C6":      impurity = 'C$^{6+}$'
    if name=="W34":     impurity = 'W$^{34+}$' 

    # Return the impurity name 
    return impurity