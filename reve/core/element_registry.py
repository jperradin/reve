from ..data import chemical_symbols, atomic_masses, atomic_names, correlation_lengths

class ElementRegistry:
    _elements = {}
    _initialized = False

    @classmethod
    def register_element(cls, element, properties):
        cls._elements[element] = properties

    @classmethod
    def register_element_component(cls, element, component):
        if element in cls._elements:
            cls._elements[element]["components"].append(component)

    @classmethod
    def get_properties(cls, element):
        if not cls._initialized:
            cls._initialize_elements()
        return cls._elements.get(element, {})

    @classmethod
    def supported_elements(cls):
        if not cls._initialized:
            cls._initialize_elements()
        return list(cls._elements.keys())
        
    @classmethod
    def _initialize_elements(cls):
        """Initialize elements from reve.data.__init__.py"""
        if cls._initialized:
            return
            
        # Register all elements from the chemical_symbols array
        for i, symbol in enumerate(chemical_symbols):
            if symbol:  # Skip empty symbol at index 0
                properties = {
                    "atomic_number": i,
                    "atomic_mass": atomic_masses[i] if i < len(atomic_masses) else 0.0,
                    "correlation_length": correlation_lengths[i] if i < len(correlation_lengths) else 0.0,
                    "chemical_symbol": symbol,
                    "atomic_name": atomic_names[i] if i < len(atomic_names) else "",
                    "components": []  # Default empty components list
                }
                cls.register_element(symbol, properties)
                
        cls._initialized = True

    