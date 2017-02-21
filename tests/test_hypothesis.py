from hypothesis import given, settings
import hypothesis.strategies as st


@given(st.integers(), st.integers())
def test_hypothesis(x, y):
    assert x + y == y + x


@settings(max_examples=200)
@given(st.lists(elements=st.integers(), min_size=2, max_size=100))
def test_google_find_pair_thing(lst):
    "google interview example that I didn't believe worked..."
    target_num = lst[0] + lst[1]
    lst.sort()
    i = 0
    j = len(lst) - 1
    found = False
    while(i < j):
        test_num = lst[i] + lst[j]
        if test_num == target_num:
            found = True
            break
        elif test_num < target_num:
            i += 1
        else:
            j -= 1
    assert found
