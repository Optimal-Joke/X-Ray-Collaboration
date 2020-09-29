# %%
import subprocess

# %%
useless_call = subprocess.run(["ls"], stdout=subprocess.PIPE, text=True)
print(useless_call.stdout)
# print(useless_call)

# %%
