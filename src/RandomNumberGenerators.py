import sys

class RANDU:
    def __init__(self, x0):
        self.x = x0
        self.a = 65539
        self.m = 2**31 - 1

    def randint(self):
        """Return a random integer in (0, 2**31)."""
        self.x = self.a * self.x % self.m
        return self.x

    def rand(self, low=0.0, high=1.0):
        """Return a random floating point number in (low, high)."""
        return low + (self.randint() / self.m) * (high - low)


class Xorshift64:
    def __init__(self, x0, a1, a2, a3):
        self.x = x0
        self.a1 = a1
        self.a2 = a2
        self.a3 = a3
        self.maxint = 2**64

    def randint_L(self):
        """Generate a random integer in (0, 2**64) using the update rule that
        starts with a left shift."""
        x = self.x ^ (self.x << self.a1) % self.maxint
        x = x ^ (x >> self.a2) % self.maxint
        self.x = (x ^ (x << self.a3)) % self.maxint
        return self.x

    def rand_L(self, low=0.0, high=1.0):
        """Generate a random floating point number in (low, high) using the update rule that
        starts with a left shift."""
        return low + (self.randint_L() / self.maxint) * (high - low)

    def randint_R(self):
        """Generate a random integer in (0, 2**64) using the update rule that
        starts with a right shift."""
        x = self.x ^ (self.x >> self.a1) % self.maxint
        x = x ^ (x << self.a2) % self.maxint
        self.x = (x ^ (x >> self.a3)) % self.maxint
        return self.x

    def rand_R(self, low=0.0, high=1.0):
        """Generate a random floating point number in (low, high) using the update rule that
        starts with a right shift."""
        return low + (self.randint_R() / self.maxint) * (high - low)


class MWC32:
    def __init__(self, x0):
        self.x = x0
        self.a = 4294957665
        self.b = 2**32
        self.m = (self.a * self.b - 2) / 2
        self.maxint = 2**64

    def randint(self):
        """Generate a random integer in [0, 2**64)."""
        self.x = (self.a * (self.x & 0xffffffff) + (self.x >> 32)) % self.maxint
        return self.x

    def rand(self, low=0.0, high=1.0):
        """Generate a random floating point number in [low, high)."""
        return low + (self.randint() / self.maxint) * (high - low)


class LCG64:
    def __init__(self, x0):
        self.x = x0
        self.a = 2862933555777941757
        self.b = 7046029254386353087
        self.m = 2**64

    def randint(self):
        """Generate a random integer in [0, 2**64)."""
        self.x = self.a * self.x % self.m
        return self.x

    def rand(self, low=0.0, high=1.0):
        """Generate a random floating point number in [low, high)."""
        return low + (self.randint() / self.m) * (high - low)


class Ran:
    def __init__(self, x0):
        """Initialize the states of the component generators C, A3R and B. The state of A1L does not
        need to be initialized, since the output (state) of C will be used as input to A1L.
        """
        self.maxint = 2**64
        self.A3R = 4101842887655102017
        self.B = 1
        self.C = x0 ^ self.A3R
        self.randint()
        self.A3R = self.C
        self.randint()
        self.B = self.A3R
        self.randint()

    def randint(self):
        """Generate a random integer in [0, 2**64)."""
        # LGC64
        self.C = (self.C * 2862933555777941757 + 7046029254386353087) % self.maxint
        # xorshift 3R
        self.A3R = self.A3R ^ (self.A3R >> 17) % self.maxint
        self.A3R = self.A3R ^ (self.A3R << 31) % self.maxint
        self.A3R = (self.A3R ^ (self.A3R >> 8)) % self.maxint
        # MWC32
        self.B = (4294957665 * (self.B & 0xffffffff) + (self.B >> 32)) % self.maxint
        # xorshift 1L, using output of LGC64 as input
        A1L = self.C ^ (self.C << 21) % self.maxint
        A1L = A1L ^ (A1L >> 35) % self.maxint
        A1L = (A1L ^ (A1L << 4)) % self.maxint
        # combine generators
        return ((A1L + self.A3R) ^ self.B) % self.maxint

    def rand(self, low=0.0, high=1.0):
        """Generate a random floating point number in [0, 1)."""
        return low + (self.randint() / self.maxint) * (high - low)


