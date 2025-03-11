import colorsys

__all__ = [
    'print_title',
    'generate_color_gradient'
]

def print_title(__version__):
    """ Print the title at REVE initialization """
    TITLE = fr'''
                                                 
`7MM"""Mq.  `7MM"""YMM `7MMF'   `7MF'`7MM"""YMM  
  MM   `MM.   MM    `7   `MA     ,V    MM    `7  
  MM   ,M9    MM   d      VM:   ,V     MM   d    
  MMmmdM9     MMmmMM       MM.  M'     MMmmMM    
  MM  YM.     MM   Y  ,    `MM A'      MM   Y  , 
  MM   `Mb.   MM     ,M     :MM;       MM     ,M 
.JMML. .JMM..JMMmmmmMMM      VF      .JMMmmmmMMM 

    '''
    print(TITLE)
    print(f"version \u279c\t {__version__}\n")


def generate_color_gradient(num_iterations):
    """ Generate a gradient of colors to update at each tqdm iteration """

    # Define the start and end colors in RGB
    start_color = (255, 0, 0)  # Red
    end_color = (0, 0, 255)    # Blue
    
    # Check if num_iterations is 0
    if num_iterations == 0:
        return [start_color]  # Return a list with only the start color

    # Check if num_iterations is 1
    elif num_iterations == 1:
        return [start_color, end_color]  # Return a list containing both start and end colors
    else:
        num_iterations += 1 
            
    # Convert RGB to HSV
    start_hsv = colorsys.rgb_to_hsv(*[x / 255.0 for x in start_color])
    end_hsv = colorsys.rgb_to_hsv(*[x / 255.0 for x in end_color])

    # Interpolate between the start and end colors
    color_gradient = []
    for i in range(num_iterations):
        ratio = i / (num_iterations - 1)
        hsv = (
            start_hsv[0] + ratio * (end_hsv[0] - start_hsv[0]),
            start_hsv[1] + ratio * (end_hsv[1] - start_hsv[1]),
            start_hsv[2] + ratio * (end_hsv[2] - start_hsv[2])
        )
        rgb = tuple(int(x * 255) for x in colorsys.hsv_to_rgb(*hsv))
        color_gradient.append(rgb)

    return color_gradient
