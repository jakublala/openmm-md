import os


def main(
        filename=None
):
    if filename is None:
        raise ValueError('Filename is required')
    
    # print("-----Running fixer.py-----")
    # os.system(f'python fixer.py --filename {filename}')
    
    # print("-----Running minimize.py-----")
    # os.system(f'python minimize.py --filename {filename}')

    print("-----Running stability.py-----")
    os.system(f'python stability.py --filename {filename}')


import fire
if __name__ == '__main__':
    fire.Fire(main)