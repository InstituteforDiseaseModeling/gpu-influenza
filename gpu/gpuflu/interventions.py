#! /usr/bin/env python3

import unittest

_MAX_DURATION=65535  # max uint16_t


class Intervention(object):

    ATTRIBUTE = 1
    NEIGHBOR = 2
    NODEMATCHING = 3
    NEIGHBORMATCHING = 4


    def __init__(self, kind, tag="unnamed", trigger="TIME", trigger_value=0, duration=_MAX_DURATION, factor=1.0, attribute=None, attribute_value=None, comment=None):
        self._type = kind
        self._tag = tag
        self._trigger = trigger.upper()
        self._trigger_value = trigger_value
        self._duration = duration
        self._factor = factor
        self._attribute = attribute if attribute else "n/a"
        self._attribute_value = attribute_value if attribute_value else 0
        self._comment = str(comment) if comment else None

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

    @property
    def comment(self):
        return self._comment

    def __str__(self):
        prefix = f"# {self.comment}\n" if self.comment else ""
        return f"{prefix}{self.tag},{self.type},{self.trigger},{self.trigger_value},{self.duration},{self.factor},{self.attribute},{self.attribute_value}"


class AttributeIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor, attribute, attribute_value, comment=None):
        super().__init__(Intervention.ATTRIBUTE, tag, trigger, trigger_value, duration, factor, attribute, attribute_value, comment)

        return


class NeighborIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor, comment=None):
        super().__init__(Intervention.NEIGHBOR, tag, trigger, trigger_value, duration, factor, comment=comment)

        return


class NodeMatchingIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor, attribute, comment=None):
        super().__init__(Intervention.NODEMATCHING, tag, trigger, trigger_value, duration, factor, attribute, comment=comment)

        return


class NeighborMatchingIntervention(Intervention):

    def __init__(self, tag, trigger, trigger_value, duration, factor, attribute, attribute_value, comment=None):
        super().__init__(Intervention.NEIGHBORMATCHING, tag, trigger, trigger_value, duration, factor, attribute, attribute_value, comment)

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
                handle.write(f"{intervention}\n")

        return


class CampaignTest(unittest.TestCase):

    def test_campaign(self):

        c = Campaign()
        c += AttributeIntervention("attribute_iv", "COUNT", 20, 365, 0.0625, "household", 0, "if count >= 20, reduce susceptibility of people in primary location 0")
        c += AttributeIntervention("attribute_iv", "COUNT", 20, 365, 0.0625, "household", 1)
        c += NeighborIntervention("neighbor_iv", "COUNT", 100, 14, 0.000030517578125, comment="if count >= 100, quarantine (2^-15) all currently infectious individuals for 14 time steps")
        c += NodeMatchingIntervention("nodematch_iv", "TIME", 60, 31, 0.125, "household")
        c += NeighborMatchingIntervention("neighbormatch_iv", "COUNT", 128, 60, 0.125, "school", 13, "if count >= 128, quarantine all households of infectious individuals")
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