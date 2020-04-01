#! /usr/bin/env python3

import unittest

_MAX_DURATION=65535  # max uint16_t


class Intervention(object):

    ATTRIBUTE = 1
    NEIGHBOR = 2
    NODEMATCHING = 3
    NEIGHBORMATCHING = 4


    def __init__(self, kind, tag="unnamed", trigger="TIME", trigger_value=0, duration=_MAX_DURATION, factor=1.0, attribute=None, attribute_value=None):
        self._type = kind
        self._tag = tag
        self._trigger = trigger.upper()
        self._trigger_value = trigger_value
        self._duration = duration
        self._factor = factor
        self._attribute = attribute if attribute else "n/a"
        self._attribute_value = attribute_value if attribute_value else 0

        self._validate()

        return

    def _validate(self):
        if self.trigger not in { "TIME", "COUNT" }:
            raise RuntimeError(f"trigger must be 'TIME' or 'COUNT', got '{self.trigger}'")
        if self.trigger_value < 0:
            raise RuntimeError(f"trigger value must be >= 0, got '{self.trigger_value}'")
        if self.duration < 1:
            raise RuntimeError(f"duration must be > 0, got '{self.duration}'")
        if self.factor < 0:
            raise RuntimeError(f"parameter modification factor must be >= 0, got '{self.factor}'")

        return

    @property
    def type(self):
        return self._type

    @property
    def tag(self):
        return self._tag

    @property
    def trigger(self):
        return self._trigger

    @property
    def trigger_value(self):
        return self._trigger_value

    @property
    def duration(self):
        return self._duration

    @property
    def factor(self):
        return self._factor

    @property
    def attribute(self):
        return self._attribute

    @property
    def attribute_value(self):
        return self._attribute_value

    def __str__(self):
        return f"{self.tag},{self.type},{self.trigger},{self.trigger_value},{self.duration},{self.factor},{self.attribute},{self.attribute_value}"


class AttributeIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor, attribute, attribute_value):
        super().__init__(Intervention.ATTRIBUTE, tag, trigger, trigger_value, duration, factor, attribute, attribute_value)

        return


class NeighborIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor):
        super().__init__(Intervention.NEIGHBOR, tag, trigger, trigger_value, duration, factor)

        return


class NodeMatchingIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor, attribute):
        super().__init__(Intervention.NODEMATCHING, tag, trigger, trigger_value, duration, factor, attribute)

        return


class NeighborMatchingIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor, attribute, attribute_value):
        super().__init__(Intervention.NEIGHBORMATCHING, tag, trigger, trigger_value, duration, factor, attribute, attribute_value)

        return


class Campaign(object):

    def __init__(self):
        self._interventions = []

        return

    @property
    def interventions(self):
        return self._interventions

    def __iadd__(self, other):
        self.interventions.append(other)
        return self

    def from_file(self, filename):

        return

    def write_file(self, filename):

        with open(filename, "wt") as handle:
            handle.write("tag,type,trigger,trigger_value,duration,factor,attribute,attribute_value"+"\n")
            for intervention in self.interventions:
                handle.write(str(intervention)+"\n")

        return


class CampaignTest(unittest.TestCase):

    def test_campaign(self):

        c = Campaign()
        a = AttributeIntervention("attribute_iv", "TIME", 30, 30, 0.25, "household", 2020)
        c += a
        n = NeighborIntervention("neighbor_iv", "COUNT", 100, 14, 0.125)
        c += n
        node = NodeMatchingIntervention("nodematch_iv", "TIME", 60, 31, 0.125, "household")
        c += node
        neighbor_match = NeighborMatchingIntervention("neighbormatch_iv", "COUNT", 128, 60, 0.125, "school", 13)
        c += neighbor_match
        c.write_file("campaign.csv")

        return

    def test_invalid_trigger_type(self):

        fn = lambda : Intervention(Intervention.ATTRIBUTE, "bad-trigger-type", "deadbeef", 100, 30, 0.5, "attribute", 2020)
        self.assertRaises(RuntimeError, fn)

        return

    def test_invalid_trigger_value(self):

        fn = lambda : Intervention(Intervention.ATTRIBUTE, "bad-trigger-type", "TIME", -2, 30, 0.5, "attribute", 2020)
        self.assertRaises(RuntimeError, fn)

        return

    def test_invalid_duration(self):

        fn = lambda : Intervention(Intervention.ATTRIBUTE, "bad-trigger-type", "TIME", 100, -30, 0.5, "attribute", 2020)
        self.assertRaises(RuntimeError, fn)

        return

    def test_invalid_factor(self):

        fn = lambda : Intervention(Intervention.ATTRIBUTE, "bad-trigger-type", "TIME", 100, 30, -0.5, "attribute", 2020)
        self.assertRaises(RuntimeError, fn)

        return



if __name__ == "__main__":
    unittest.main()