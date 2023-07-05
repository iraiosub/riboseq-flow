import dnaio
import gzip
import argparse


def break_to_list(sequence, length):
	list = []
	i = -1
	while True:
		i += 1
		if i + length > len(sequence):
			break
		list.append(sequence[i:i+length])

	return list


def main():

	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required=True)
	parser.add_argument("-o", "--output", required=True)
	parser.add_argument("-r", "--reference", required=True)
	parser.add_argument("-l", "--length", type=int, default=14)
	args = parser.parse_args()

	rRNA_list = []

	with open(args.reference) as fasta:
		for i, line in enumerate(fasta):
			if i % 2 == 1:
				rRNA_list += break_to_list(line.rstrip(), args.length)

	rRNA_set = set(rRNA_list)

	print(rRNA_set)

	out_list = []

	n = 0
	good = 0
	with dnaio.open(args.input) as file, gzip.open(args.output, 'wb') as out:
		for record in file:
			n += 1
			s = str(record.sequence)
			print(s)
			seq_list = break_to_list(s, args.length)

			union_n = set(seq_list).intersection(rRNA_set)

			print(seq_list)

			print(union_n)
			print(len(union_n))

			if n % 100 == 0:
				print(n)
				break

			if len(union_n) == 0:
				good += 1
				out_list += [str(record.name), str(record.sequence), "+", str(record.qualities)]

		print(str(good) + " passed filter out of total = " + str(n))

		out.write(('\n'.join(out_list)).encode())


if __name__ == "__main__":
	main()