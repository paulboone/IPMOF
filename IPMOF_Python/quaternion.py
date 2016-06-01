class Quaternion(object):
    def __init__(self, input):
        if isinstance(input, list):
            self.w = input[0]
            self.x = input[1]
            self.y = input[2]
            self.z = input[3]
            self.q = [self.x, self.y, self.z]
        else:
            raise TypeError('Input type not supported. Use list [w, x, y, z]')

    def xyz(self):
        return [self.x, self.y, self.z]

    def __mul__(self, quat2):
        """
        Multiply quaternion by another.

        Example usage::

          >>> q1 = Quat((20,30,40))
          >>> q2 = Quat((30,40,50))
          >>> (q1 * q2).equatorial
          array([ 349.73395729,   76.25393056,  127.61636727])

        :returns: product q1 * q2
        :rtype: Quat

        """

        q1 = self
        q2 = quat2

        w3 = q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
        x3 = q1.x * q2.w + q1.w * q2.x - q1.z * q2.y + q1.y * q2.z
        y3 = q1.y * q2.w + q1.z * q2.x + q1.w * q2.y - q1.x * q2.z
        z3 = q1.z * q2.w - q1.y * q2.x + q1.x * q2.y + q1.w * q2.z

        return Quaternion([w3, x3, y3, z3])

    def __div__(self, quat2):
        """
        Divide one quaternion by another.

        Example usage::

         >>> q1 = Quat((20,30,40))
         >>> q2 = Quat((30,40,50))
         >>> q = q1 / q2

        Performs the operation as q1 * inverse q2

        :returns: product q1 * inverse q2
        :rtype: Quat

        """
        return self * quat2.inv()

    def inv(self):
        """
        Invert the quaternion

        :returns: inverted quaternion
        :rtype: Quat
        """
        norm = self.w ** 2 + self.x ** 2 + self.y ** 2 + self.z ** 2

        return Quaternion([self.w / norm, -self.x / norm, -self.y / norm, -self.z / norm])
