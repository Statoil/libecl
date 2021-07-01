#  Copyright (C) 2014.  Equinor ASA, Norway.
#
#  This file is part of ERT - Ensemble based Reservoir Tool.
#
#  ERT is free software: you can redistribute it and/or modify it under the
#  terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later
#  version.
#
#  ERT is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE.
#
#  See the GNU General Public License at <http://www.gnu.org/licenses/gpl.html>
#  for more details.

from __future__ import print_function


class FaultSegment(object):

    def __init__(self, C1, C2):
        self.__C1 = C1
        self.__C2 = C2
        self.__next_segment = None


    def __eq__(self, other):
        s = self.c1, self.c2
        o = other.c1, other.c2
        o_flipped = other.c2, other.c1
        return s == o or s == o_flipped

    def __hash__(self):
        return hash(hash(self.__C1) + hash(self.__C2) + hash(self.__next_segment))

    def get_corners(self):
        return (self.__C1, self.__C2)

    def joins(self, other):
        if self.__C1 == other.__C1:
            return True
        if self.__C1 == other.__C2:
            return True
        if self.__C2 == other.__C1:
            return True
        if self.__C2 == other.__C2:
            return True

        return False

    def get_c1(self):
        return self.__C1

    @property
    def c1(self):
        return self.__C1

    def get_c2(self):
        return self.__C2

    @property
    def c2(self):
        return self.__C2

    def swap(self):
        C1 = self.__C1
        self.__C1 = self.__C2
        self.__C2 = C1


    def __repr__(self):
        return "%d -> %d" % (self.__C1, self.__C2)



class SegmentMap(object):

    def __init__(self):
        self.__segment_map = {}
        self.__count_map = {}

    def __len__(self):
        return len(self.__segment_map)

    def __str__(self):
        return self.__segment_map.__str__()

    def verify(self):
        for (C, count) in self.__count_map.iteritems():
            if count > 0:
                d = self.__segment_map[C]
                if len(d) != count:
                    print('CornerPoint:%d  count:%d  len(d):%d map:%s' % (C, count, len(d), d))
                assert len(d) == count
            else:
                assert self.__segment_map.get(C) is None


    def add_segment(self, segment):
        (C1,C2) = segment.get_corners()
        if C1 not in self.__segment_map:
            self.__segment_map[C1] = {}
            self.__count_map[C1] = 0
        if C2 not in self.__segment_map:
            self.__segment_map[C2] = {}
            self.__count_map[C2] = 0

        if C2 not in self.__segment_map[C1]:
            self.__segment_map[C1][C2] = segment
            self.__count_map[C1] += 1

        if C1 not in self.__segment_map[C2]:
            self.__segment_map[C2][C1] = segment
            self.__count_map[C2] += 1



    def del_segment(self, segment):
        (C1,C2) = segment.get_corners()
        self.__count_map[C1] -= 1
        self.__count_map[C2] -= 1
        del self.__segment_map[C1][C2]
        del self.__segment_map[C2][C1]

        if len(self.__segment_map[C1]) == 0:
            del self.__segment_map[C1]

        if len(self.__segment_map[C2]) == 0:
            del self.__segment_map[C2]


    def pop_start(self):
        end_segments = []
        for (C, count) in self.__count_map.items():
            if count == 1:
                end_segments.append(list(self.__segment_map[C].values())[0])

        start_segment = end_segments[0]
        self.del_segment(start_segment)
        return start_segment

    def pop_next(self, segment):
        (C1,C2) = segment.get_corners()
        if self.__count_map[C1] >= 1:
            next_segment = list(self.__segment_map[C1].values())[0]
        elif self.__count_map[C2] >= 1:
            next_segment = list(self.__segment_map[C2].values())[0]
        else:
            next_segment = None

        if next_segment:
            self.del_segment(next_segment)
        return next_segment


    def print_content(self):
        for d in self.__segment_map.values():
            for (C,S) in d.iteritems():
                print(S)
