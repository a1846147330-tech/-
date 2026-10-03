"""Write a deterministic synthetic OLE/BIFF XLS using only Python stdlib."""
from pathlib import Path
import math
import random
import struct


def record(tag, payload=b''):
    return struct.pack('<HH', tag, len(payload)) + payload


def generate(path):
    rng = random.Random(1)
    bof = lambda kind: record(0x809, struct.pack('<HHHHII', 0x600, kind, 0x0DBB, 1997, 0x41, 6))
    sheet = bof(0x10)
    sheet += record(0x200, struct.pack('<IIHHH', 0, 181, 0, 9, 0))
    levels = [1.0] * 8
    for row in range(1, 181):
        common = rng.gauss(0, 0.001)
        values = [46000 + row]
        for j in range(8):
            levels[j] += 0.00015 + common + rng.gauss(0, 0.001 + j * 0.00008) + 0.0002 * math.sin(row / (4 + j))
            values.append(levels[j])
        for col, value in enumerate(values):
            sheet += record(0x203, struct.pack('<HHHd', row, col, 0, value))
    sheet += record(0x0A)
    prefix = bof(5) + record(0x22, struct.pack('<H', 0))
    bound = record(0x85, struct.pack('<IBBBB', 0, 0, 0, 9, 0) + b'Synthetic')
    offset = len(prefix) + len(bound) + 4
    bound = record(0x85, struct.pack('<IBBBB', offset, 0, 0, 9, 0) + b'Synthetic')
    book = prefix + bound + record(0x0A) + sheet
    sectors = (len(book) + 511) // 512
    directory_id, fat_id = sectors, sectors + 1
    if fat_id >= 128:
        raise ValueError('Example exceeds single FAT sector')
    header = bytearray(512)
    header[:8] = bytes.fromhex('D0CF11E0A1B11AE1')
    struct.pack_into('<HHHHH', header, 24, 0x3e, 3, 0xfffe, 9, 6)
    struct.pack_into('<IIIIIIIII', header, 40, 0, 1, directory_id, 0, 4096, 0xfffffffe, 0, 0xfffffffe, 0)
    struct.pack_into('<109I', header, 76, fat_id, *([0xffffffff] * 108))

    def entry(name, kind, child=0xffffffff, start=0xfffffffe, size=0):
        item = bytearray(128)
        encoded = (name + '\0').encode('utf-16le')
        item[:len(encoded)] = encoded
        struct.pack_into('<HBBIII', item, 64, len(encoded), kind, 1, 0xffffffff, 0xffffffff, child)
        struct.pack_into('<IQ', item, 116, start, size)
        return item

    directory = entry('Root Entry', 5, child=1) + entry('Workbook', 2, start=0, size=len(book)) + bytes(256)
    fat = [0xffffffff] * 128
    for i in range(sectors):
        fat[i] = i + 1 if i + 1 < sectors else 0xfffffffe
    fat[directory_id], fat[fat_id] = 0xfffffffe, 0xfffffffd
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + book.ljust(sectors * 512, b'\0') + directory + struct.pack('<128I', *fat))
    print(path)


if __name__ == '__main__':
    generate(Path(__file__).resolve().parents[1] / 'examples/synthetic.xls')
