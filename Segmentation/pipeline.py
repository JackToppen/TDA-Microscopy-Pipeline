#!pip install pandas 
import os
from google.colab import drive
drive.mount('/content/drive', force_remount=True)

input_folder = "/content/drive/My Drive/segmentation/seleipiri"
output_folder = "/content/drive/My Drive/segmentation/seleipiri/output"
csv_output = "intensities_bgc.csv"
csv_output_path = os.path.join(output_folder,csv_output)

dapi = "D28_So10C3_03_DAPI.tiff"
dapi_path = os.path.join(input_folder,dapi)
ch_paths = [os.path.join(input_folder,"D28_So10C3_03_SOX2.tiff"),
            os.path.join(input_folder,"D28_So10C3_03_KI67.tiff"),
            os.path.join(input_folder,"D28_So10C3_03_TBR1.tiff")]

num_sizes = 10
min_rad_nucleus = 5 # approximate radius of the smallest cell in pixels
max_rad_nucleus = 25 # approximate radius of the largest cell in pixels
dapi_threshold = 500 # filter out dim nuclei (for 16-bit image max is 65536, so 500 should be low enough)
shape_channel_num = 1 # put 0 to use dapi for morphology metrics calculations

dpi = 600 # output image quality
# Potential error: 'NoneType' object has no attribute 'shape'. It means the file was not found
#Remember to change the name of the output path

image = Image(dapi_path,16,ch_paths)
image.read()
image.segment_pearson(min_rad_nucleus,max_rad_nucleus,num_sizes,channel=shape_channel_num)
image.getcellprops(output_folder,dpi)
image.overlay_cells()
print("overlaid mask and image: blue is leftover dapi staining (false negative?),")
print("green is false positive cells, light blue is correct segmentation")
mask = np.zeros((image.h,image.w,3))
mask[:,:,1] = np.array(image.mask,dtype="uint16")
mask[:,:,2] = image.dapi/np.mean(image.dapi)
fig=plt.figure(figsize=(10,10))
plt.imshow(mask)
plt.show()
image.df.to_csv(csv_output_path)
save_df = image.df
print("csv file with cell intensities is saved to",csv_output_path)
