import click
import sys
import os
import shutil
import json

__cdir__ = os.path.dirname(os.path.abspath(__file__))

@click.command(context_settings=dict(
    ignore_unknown_options=True,
    allow_extra_args=True
))
@click.pass_context
@click.option('--brick', '-b', help='Brick name')
@click.option('--force', is_flag=True, help='Remove existing brick with the same name')
def create_brick(ctx, brick=None, force=False):
    
    brick = brick.lower().replace(" ", "_").replace("-", "_")
    
    if not brick:
        print("Cannot create the brick. The brick name is required [OPTION --brick]")
        return
    
    skeleton_dir = os.path.join(__cdir__, "../../skeleton")
    dest_dir = os.path.join(__cdir__, "../../", brick)

    if os.path.exists(dest_dir):
        if not force:
            print("Cannot create the brick. A brick already exists with this name [OPTION --force allows overriding any existing brick with the same name.]")
            return
        else:
            shutil.rmtree(dest_dir)
    
    shutil.copytree(
        skeleton_dir, 
        dest_dir,
        dirs_exist_ok=True
    )
    
    if os.path.exists(os.path.join(dest_dir, brick)):
        shutil.rmtree(dest_dir)
        print(f"The brick name '{brick}' is not valid")
        return
        
    shutil.move(
        os.path.join(dest_dir, "skeleton"), 
        os.path.join(dest_dir, brick)
    )

    # remove .git folder
    shutil.rmtree(os.path.join(dest_dir, ".git"))

    
    # update settings.json
    settings_file = os.path.join(dest_dir, "settings.json")
    with open(settings_file, 'r') as f:
        settings = json.load(f)
        settings["name"]                = brick

    with open(settings_file, 'w') as f:
        json.dump(settings, f, indent=4)

    #replace all words 'skeleton' in settings.json
    with open(settings_file, 'r') as f:
        text = f.read()
        text = text.replace("skeleton", brick)

    with open(settings_file, 'w') as f:
        f.write(text)

    #replace all words 'skeleton' in app.py
    file = os.path.join(dest_dir, brick, "./app.py")
    with open(file, 'r') as f:
        text = f.read()
        text = text.replace("skeleton", brick)
    with open(file, 'w') as f:
        f.write(text)
        
    #replace all words 'skeleton' in README.md
    file = os.path.join(dest_dir, "./README.md")
    with open(file, 'r') as f:
        text = f.read()
        text = text.replace("skeleton", brick)
        text = text.replace("Skeleton", brick.title())
    with open(file, 'w') as f:
        f.write(text)
    
    print("The brick was successfully created.")

if __name__ == "__main__":
    create_brick()