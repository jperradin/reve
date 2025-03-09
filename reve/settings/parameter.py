class Setting:
    def __init__(self, name, value):
        self.name = name
        self.value = value
    
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, value):
        self._name = value
    
    def set_value(self, new_value):
        self.value = new_value
    
    def __str__(self):
        return f"{self.name}: {self.value}"
    
    def __repr__(self):
        return f"{self.name}: {self.value}"