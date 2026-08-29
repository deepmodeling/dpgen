import unittest

from dpgen.data.surf import _fixed_layers_enabled


class TestFixedLayerOptions(unittest.TestCase):
    def test_both_options_enable_fixed_layers(self):
        self.assertTrue(
            _fixed_layers_enabled({"fix_layers": 2, "total_layers": 5})
        )

    def test_no_options_disable_fixed_layers(self):
        self.assertFalse(_fixed_layers_enabled({}))

    def test_partial_options_are_rejected(self):
        for settings in ({"fix_layers": 2}, {"total_layers": 5}):
            with self.subTest(settings=settings):
                with self.assertRaisesRegex(ValueError, "provided together"):
                    _fixed_layers_enabled(settings)


if __name__ == "__main__":
    unittest.main()
