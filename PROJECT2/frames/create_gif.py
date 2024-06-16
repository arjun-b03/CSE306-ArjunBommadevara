from PIL import Image
import os
import re


image_folder = '/Users/arjunbommadevara/Desktop/CSE306/project1/CSE306-Project-1/PROJECT2/frames'

output_gif = 'animation.gif'



filenames = os.listdir(image_folder)

def extract_number(filename):
    match = re.search(r'(\d+)', filename)
    return int(match.group(1)) if match else 0

sorted_filenames = sorted(filenames, key=extract_number)

print(sorted_filenames)

image_files = [os.path.join(image_folder, f) for f in sorted_filenames if f.startswith('frame_') and f.endswith('.png')]


images = [Image.open(image) for image in image_files]


images[0].save(output_gif, save_all=True, append_images=images[1:], duration=50, loop=0)
