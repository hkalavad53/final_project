import sys

def calculate_lowercase_percentage(text):
    total_characters = len(text)
    lowercase_count = sum(1 for char in text if char.islower())
    lowercase_percentage = (lowercase_count / total_characters) * 100
    return lowercase_percentage

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <file_path>")
        return
    
    file_path = sys.argv[1]
    
    try:
        with open(file_path, 'r') as file:
            content = file.read()
            lines = content.split('\n')
            
            total_lowercase_count = 0
            total_characters_count = 0
            
            for line in lines:
                if not line.startswith(">"):
                    total_lowercase_count += sum(1 for char in line if char.islower())
                    total_characters_count += len(line)
            
            lowercase_percentage = (total_lowercase_count / total_characters_count) * 100
            print(f"Percentage of lowercase letters: {lowercase_percentage:.2f}%")
    except FileNotFoundError:
        print("File not found.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()

