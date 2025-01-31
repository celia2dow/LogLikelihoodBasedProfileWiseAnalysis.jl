function get_darker_colour(colour_input, amount=0.2)
    # Takes a Symbol or RGB colour in as input and returns a slightly (by 20%) darker colour
    
    # Dowling (2024)

    # Convert the colour to RGB if it is a symbol
    if colour_input isa Symbol
        colour_rgb = parse(Colorant, colour_input)
    else
        colour_rgb = colour_input
    end
    
    # Convert the colour to HSV, reduce the value (brightness), then convert back to RGB
    hsv = convert(HSV, colour_rgb)
    darker_hsv = HSV(hsv.h, hsv.s, max(0, hsv.v - amount))
    colour_input_darker = convert(RGB, darker_hsv)
    
    return colour_input_darker
end