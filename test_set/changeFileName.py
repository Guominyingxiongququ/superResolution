import sys
import getopt
import os

def main():
    num = 0 
    for filename in os.listdir("."):
        if filename.startswith("COCO_test2014_"):
            num = num + 1
            os.rename(filename, "img_" + str(num) + ".jpg")

if __name__ == "__main__":
    main()