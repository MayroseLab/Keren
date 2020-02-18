import os, argparse

if __name__ == '__main__':

    # process input from command line
    parser = argparse.ArgumentParser(description='removes spaces from all the files in a directory')
    parser.add_argument('--input_dir', '-i', help='directory of files from which spaces should be removed',
                        required=True)

    args = parser.parse_args()
    input_dir = args.input_dir

    # remove spaces from simulated MSAs (because it makes Bio++ shout)
    print("removing spaces from taxa names in the simulated MSAs")
    for file in os.listdir(input_dir):
        with open(input_dir+file, "r") as input:
            file_content = input.read()
            file_content = file_content.replace(" ", "")
        with open(input_dir+file, "w") as output:
            output.write(file_content)